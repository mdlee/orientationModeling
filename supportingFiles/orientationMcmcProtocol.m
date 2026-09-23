function [chains, stats, diagnostics, info] = orientationMcmcProtocol(cfg)
%ORIENTATIONMCMCPROTOCOL Chunked JAGS sampling without chain surgery.
%
% Required cfg fields:
%   modelFile, modelName, data, params, nStimuli, stimuli, makeInits,
%   lastStateInits, sizeField
%
% makeInits is () -> struct. lastStateInits(chains) -> 1 x nChains cell.
% sizeField is a coda field used for sample/chain counts ('sigma' or 'sigma_1').
%
% Optional: modules, nChains, nBurninAdapt, nBurninRestart, nThin,
% nDrawsWanted, nSamplesChunk, nSamplesMax, rHatMax, essMin, doParallel,
% rngSeed, piVal, restartJitter, extraDiagPrefixes.

cfg = fillDefaults(cfg);

if isfield(cfg, 'chains') && ~isempty(cfg.chains)
   chains = dropWrapMus(cfg.chains, cfg.nStimuli);
   [chains, nFlips] = alignMuWraps(chains, cfg.stimuli, cfg.piVal, ...
      cfg.nStimuli, cfg.sizeField);
   if nFlips > 0
      fprintf('Wrap-aligned %d additional chain-by-bin mu shifts after load\n', nFlips);
   end
   diagnostics = diagnosePosterior(chains, cfg);
   printDiagnostics(diagnostics);
   stats = [];
   info = [];
   if isfield(cfg, 'stats'), stats = cfg.stats; end
   if isfield(cfg, 'info'), info = cfg.info; end
   return
end

rng(cfg.rngSeed);

chains = [];
stats = [];
info = [];
chunk = 0;
nBurnin = cfg.nBurninAdapt;
nTake = ceil(cfg.nDrawsWanted / cfg.nChains);
initArg = cfg.makeInits;

while true
   chunk = chunk + 1;
   fprintf('MCMC chunk %d: %d chains, burn-in %d, kept %d, thin %d\n', ...
      chunk, cfg.nChains, nBurnin, nTake, cfg.nThin);

   [statsChunk, chainsChunk, ~, infoChunk] = callJags(cfg, initArg, ...
      nBurnin, nTake, cfg.nChains, cfg.modelName, chunk, 0);

   if isempty(chains)
      chains = chainsChunk;
      stats = statsChunk;
      info = infoChunk;
   else
      chains = appendChains(chains, chainsChunk);
      stats = statsChunk;
      info = infoChunk;
   end

   chains = dropWrapMus(chains, cfg.nStimuli);
   [chains, nFlips] = alignMuWraps(chains, cfg.stimuli, cfg.piVal, ...
      cfg.nStimuli, cfg.sizeField);
   fprintf('Wrap-aligned %d chain-by-bin mu shifts\n', nFlips);

   chains = restartLostChains(chains, cfg, chunk);

   diagnostics = diagnosePosterior(chains, cfg);
   printDiagnostics(diagnostics);
   nKept = size(chains.(cfg.sizeField), 1);

   if diagnostics.converged
      fprintf('Stop rule passed with %d draws/chain (%d joint).\n', ...
         nKept, nKept * cfg.nChains);
      break
   end
   if nKept >= cfg.nSamplesMax
      warning('orientationMcmcProtocol:noConverge', ...
         ['Reached per-chain cap (%d) without passing the stop rule. ', ...
         'Saving all %d chains anyway; inspect diagnostics.'], ...
         cfg.nSamplesMax, cfg.nChains);
      break
   end

   nTake = min(cfg.nSamplesChunk, cfg.nSamplesMax - nKept);
   nBurnin = cfg.nBurninRestart;
   setChainInitPool(cfg.lastStateInits(chains));
   initArg = @nextChainInit;
end

disp('Convergence statistics (Rhat >= 1.05):')
grtable(chains, cfg.rHatMax)
disp('Descriptive statistics for all chains:')
codatable(chains);
diagnostics = diagnosePosterior(chains, cfg);
end

function cfg = fillDefaults(cfg)
req = {'modelFile', 'modelName', 'data', 'params', 'nStimuli', ...
   'stimuli', 'makeInits', 'lastStateInits', 'sizeField'};
for i = 1:numel(req)
   if ~isfield(cfg, req{i})
      error('orientationMcmcProtocol:missing', 'cfg.%s is required', req{i});
   end
end
if ~isfield(cfg, 'modules'), cfg.modules = {}; end
if ~isfield(cfg, 'nChains'), cfg.nChains = 8; end
if ~isfield(cfg, 'nBurninAdapt'), cfg.nBurninAdapt = 7000; end
if ~isfield(cfg, 'nBurninRestart'), cfg.nBurninRestart = 2000; end
if ~isfield(cfg, 'nThin'), cfg.nThin = 10; end
if ~isfield(cfg, 'nDrawsWanted'), cfg.nDrawsWanted = 1e4; end
if ~isfield(cfg, 'nSamplesChunk'), cfg.nSamplesChunk = 250; end
if ~isfield(cfg, 'nSamplesMax'), cfg.nSamplesMax = 2500; end
if ~isfield(cfg, 'rHatMax'), cfg.rHatMax = 1.05; end
if ~isfield(cfg, 'essMin'), cfg.essMin = 400; end
if ~isfield(cfg, 'doParallel'), cfg.doParallel = 1; end
if ~isfield(cfg, 'rngSeed'), cfg.rngSeed = 1; end
if ~isfield(cfg, 'piVal'), cfg.piVal = 3.1415; end
if ~isfield(cfg, 'restartJitter'), cfg.restartJitter = @(t) t; end
if ~isfield(cfg, 'extraDiagPrefixes'), cfg.extraDiagPrefixes = {}; end
end

function [stats, chains, diagnostics, info] = callJags(cfg, initArg, ...
   nBurnin, nSamples, nChains, tag, chunk, chainId)
logTag = tag;
if chainId > 0
   logTag = sprintf('%s_r%d', tag, chainId);
end
args = { ...
   'model'           , cfg.modelFile                              , ...
   'data'            , cfg.data                                   , ...
   'outputname'      , 'samples'                                  , ...
   'init'            , initArg                                    , ...
   'datafilename'    , logTag                                     , ...
   'initfilename'    , logTag                                     , ...
   'scriptfilename'  , logTag                                     , ...
   'logfilename'     , sprintf('/tmp/%s_c%d', logTag, chunk)      , ...
   'nchains'         , nChains                                    , ...
   'nburnin'         , nBurnin                                    , ...
   'nsamples'        , nSamples                                   , ...
   'monitorparams'   , cfg.params                                 , ...
   'thin'            , cfg.nThin                                  , ...
   'workingdir'      , sprintf('/tmp/%s_c%d', logTag, chunk)      , ...
   'verbosity'       , 0                                          , ...
   'saveoutput'      , true                                       , ...
   'allowunderscores', 1                                          , ...
   'parallel'        , cfg.doParallel                             , ...
   'seed'            , cfg.rngSeed + 17 * chunk + chainId         };
if ~isempty(cfg.modules)
   args = [args, {'modules', cfg.modules}];
end
[stats, chains, diagnostics, info] = callbayes('jags', args{:});
end

function chains = dropWrapMus(chains, nStimuli)
fn = fieldnames(chains);
for i = 1:numel(fn)
   if strncmp(fn{i}, 'mu_', 3)
      idx = sscanf(fn{i}, 'mu_%d');
      if ~isempty(idx) && idx > nStimuli
         chains = rmfield(chains, fn{i});
      end
   end
end
end

function [chains, nFlips] = alignMuWraps(chains, stimuli, piVal, nStimuli, sizeField)
nFlips = 0;
nCh = size(chains.(sizeField), 2);
for j = 2:nStimuli
   fname = sprintf('mu_%d', j);
   if ~isfield(chains, fname)
      continue
   end
   X = chains.(fname);
   sLoc = stimuli(j);
   for c = 1:nCh
      m = circularMeanHalf(X(:, c), piVal);
      cands = [m, m - piVal, m + piVal];
      [~, k] = min(abs(cands - sLoc));
      if k == 1
         continue
      end
      nFlips = nFlips + 1;
      X(:, c) = X(:, c) + (cands(k) - m);
   end
   chains.(fname) = X;
end
end

function mu0 = circularMeanHalf(vals, piVal)
phi = 2 * vals;
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
mu0 = mod(phi0 / 2, piVal);
end

function d = halfCircleDist(a, b, piVal)
d = abs(a - b);
d = min(d, piVal - d);
end

function chains = appendChains(a, b)
chains = a;
fn = fieldnames(a);
for i = 1:numel(fn)
   if isfield(b, fn{i})
      chains.(fn{i}) = [a.(fn{i}); b.(fn{i})];
   end
end
end

function setChainInitPool(inits)
nextChainInit(inits);
end

function s = nextChainInit(inits)
persistent pool k
if nargin >= 1
   pool = inits;
   k = 0;
   s = [];
   return
end
k = k + 1;
s = pool{k};
end

function chains = restartLostChains(chains, cfg, chunk)
nCh = size(chains.(cfg.sizeField), 2);
nKept = size(chains.(cfg.sizeField), 1);
err = zeros(1, nCh);
for c = 1:nCh
   acc = 0;
   for j = 1:cfg.nStimuli
      fname = sprintf('mu_%d', j);
      if ~isfield(chains, fname)
         continue
      end
      m = circularMeanHalf(chains.(fname)(:, c), cfg.piVal);
      acc = acc + halfCircleDist(m, cfg.stimuli(j), cfg.piVal)^2;
   end
   err(c) = acc;
end
med = median(err);
lost = find(err > max(2.5 * med, med + 0.25));
if isempty(lost)
   return
end
[~, best] = min(err);
fprintf('Restarting chain(s) %s from chain %d (map error %.3f vs median %.3f)\n', ...
   mat2str(lost), best, max(err(lost)), med);
baseInits = cfg.lastStateInits(chains);
template = cfg.restartJitter(baseInits{best});
for i = 1:numel(lost)
   c = lost(i);
   [~, chainsNew, ~, ~] = callJags(cfg, @()template, ...
      cfg.nBurninRestart, nKept, 1, cfg.modelName, chunk, c);
   chainsNew = dropWrapMus(chainsNew, cfg.nStimuli);
   [chainsNew, ~] = alignMuWraps(chainsNew, cfg.stimuli, cfg.piVal, ...
      cfg.nStimuli, cfg.sizeField);
   fn = fieldnames(chains);
   for f = 1:numel(fn)
      if isfield(chainsNew, fn{f}) && size(chainsNew.(fn{f}), 2) == 1
         chains.(fn{f})(:, c) = chainsNew.(fn{f})(:, 1);
      end
   end
end
end

function diag = diagnosePosterior(chains, cfg)
names = monitoredNames(chains, cfg.nStimuli, cfg.extraDiagPrefixes);
nNames = numel(names);
rhat = nan(nNames, 1);
ess = nan(nNames, 1);
essTail = nan(nNames, 1);
for i = 1:nNames
   x = chains.(names{i});
   rhat(i) = splitRhat(x);
   try
      ess(i) = gelmanrubin(x, 0, 1, 'neff');
   catch
      ess(i) = NaN;
   end
   essTail(i) = tailEss(x);
end
nKept = size(chains.(cfg.sizeField), 1);
nCh = size(chains.(cfg.sizeField), 2);
diag.names = names;
diag.rhat = rhat;
diag.ess = ess;
diag.essTail = essTail;
diag.maxRhat = max(rhat);
diag.minEss = min(ess);
diag.minEssTail = min(essTail);
diag.nJoint = nKept * nCh;
diag.nChains = nCh;
[diag.maxRhatName, idxR] = worstName(names, rhat, @max);
[diag.minEssName, idxE] = worstName(names, ess, @min);
diag.converged = diag.maxRhat < cfg.rHatMax ...
   && diag.minEss >= cfg.essMin ...
   && diag.minEssTail >= cfg.essMin ...
   && diag.nJoint >= cfg.nDrawsWanted ...
   && nCh >= 4;
if isempty(idxR)
   diag.maxRhatName = '';
end
if isempty(idxE)
   diag.minEssName = '';
end
end

function names = monitoredNames(chains, nStimuli, extraPrefixes)
fn = fieldnames(chains);
keep = {};
for i = 1:numel(fn)
   name = fn{i};
   if strncmp(name, 'mu_', 3)
      idx = sscanf(name, 'mu_%d');
      if isempty(idx) || idx < 2 || idx > nStimuli
         continue
      end
   elseif strcmp(name, 'sigma') || strncmp(name, 'sigma_', 6)
      % keep
   elseif strcmp(name, 'psi') || strncmp(name, 'psi_', 4)
      % keep
   else
      hit = false;
      for p = 1:numel(extraPrefixes)
         if strcmp(name, extraPrefixes{p}) || strncmp(name, [extraPrefixes{p} '_'], numel(extraPrefixes{p}) + 1)
            hit = true;
            break
         end
      end
      if ~hit
         continue
      end
   end
   if std(chains.(name)(:)) < 1e-12
      continue
   end
   keep{end+1} = name; %#ok<AGROW>
end
names = keep;
end

function r = splitRhat(x)
[n, m] = size(x);
n2 = floor(n / 2);
if n2 < 50 || m < 2
   r = gelmanrubin(x, 0, 1, 'rhat');
   return
end
xsplit = [x(1:n2, :) x((n - n2 + 1):n, :)];
r = gelmanrubin(xsplit, 0, 1, 'rhat');
end

function ess = tailEss(x)
q = prctile(x(:), [5 95]);
essLo = indicatorEss(x, x <= q(1));
essHi = indicatorEss(x, x >= q(2));
ess = min(essLo, essHi);
end

function ess = indicatorEss(x, mask)
z = double(mask);
if std(z(:)) < 1e-12
   ess = size(x, 1) * size(x, 2);
   return
end
try
   ess = gelmanrubin(z, 0, 1, 'neff');
catch
   ess = NaN;
end
end

function [name, idx] = worstName(names, vals, fun)
idx = [];
name = '';
ok = isfinite(vals);
if ~any(ok)
   return
end
[~, loc] = fun(vals(ok));
idxMap = find(ok);
idx = idxMap(loc);
name = names{idx};
end

function printDiagnostics(diag)
fprintf(['Diagnostics: max split-Rhat = %.4f (%s); min bulk ESS = %.0f (%s); ', ...
   'min tail ESS = %.0f; joint draws = %d; converged = %d\n'], ...
   diag.maxRhat, diag.maxRhatName, diag.minEss, diag.minEssName, ...
   diag.minEssTail, diag.nJoint, diag.converged);
end
