%% Memory reproduction MCMC (wrap-copy Gaussian with swaps, automated)
%
% Implements the shareable sampler for the wrap-copy swap model:
%   8 chains kept in full, data-informed inits, no xi/yP monitors,
%   wrap alignment of mu, split-Rhat / ESS gates, 10k joint draws for
%   Savage-Dickey on omega. See memoryReproductionMCMC.pdf.
%
% Storage: storage/memoryReproductionMCMC_tomicBays_jags.mat
% (does not overwrite memoryReproduction_tomicBays_jags.mat)
%
% Trinity's JAGS wrapper has no separate adapt command; nBurnin = 7000 is
% ~2k JAGS adaptation + 5k post-adapt burn-in.

clear; close all;
preLoad = true;
printFigures = true;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

modelDir = './';
modelName = 'memoryReproductionMCMC';
engine = 'jags';

dataList = {...
  'tomicBaysMemory'; ...
  };

%% constants
pi = 3.1415;
load pantoneColors pantone
fontSize = 18;
CI = [2.5 97.5];

% MCMC protocol
nChains         = 8;
nBurninAdapt    = 7000;   % ~2k adapt + 5k burn-in (trinity has no adapt knob)
nBurninRestart  = 2000;   % short re-adapt when extending / restarting a chain
nThin           = 10;
nDrawsWanted    = 1e4;    % joint posterior draws for Savage-Dickey
nSamplesChunk   = 250;
nSamplesMax     = 2500;   % per chain cap
rHatMax         = 1.05;
essMin          = 400;
doParallel      = 1;
rngSeed         = 1;

nSamplesTarget = ceil(nDrawsWanted / nChains);

% loop over data
for dataIdx = 1:numel(dataList)
  dataName = dataList{dataIdx};
  switch dataName

    case 'tomicBaysMemory'
      dataDir = '../data/';
      dataName = 'tomicBays';
      load([dataDir dataName], 'dm');

      [~, ~, setSize] = unique(dm.setSize, 'stable');
      y = dm.response;
      nTrials = dm.nTrials;
      nStimuli = dm.nStimuli;
      s = [dm.tIdx dm.nIdx];
      s(isnan(s)) = 1;
      [~, maxPresented] = size(s);

      % order nontargets in s by similarity to target
      for t = 1:nTrials
        vals = dm.nontarget(t, 1:(dm.setSize(t)-1));
        [~, srt] = sort([0 min(abs(dm.target(t)-vals), pi-abs(dm.target(t)-vals))], 'ascend');
        s(t, 1:dm.setSize(t)) = s(t, srt);
      end

  end

  %% sampling from graphical model
  params = {'mu', 'sigma', 'omega3', 'omega6', 'omega3prior', 'omega6prior'};

  data = struct(...
    'setSize'  , setSize  , ...
    's'        , s        , ...
    'y'        , y        , ...
    'nStimuli' , nStimuli , ...
    'nTrials'  , nTrials  , ...
    'maxPresented', maxPresented);

  stimuli = dm.stimuli(:);
  generator = @()makeInits(stimuli, nTrials, pi);

  fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

  if preLoad && isfile(sprintf('storage/%s', fileName))
    fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
    load(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
  else
    rng(rngSeed);
    chains = [];
    stats = [];
    info = [];
    chunk = 0;
    nBurnin = nBurninAdapt;
    nTake = nSamplesTarget;
    initArg = generator;
    tic;

    while true
      chunk = chunk + 1;
      fprintf('MCMC chunk %d: %d chains, burn-in %d, kept %d, thin %d\n', ...
        chunk, nChains, nBurnin, nTake, nThin);

      [statsChunk, chainsChunk, ~, infoChunk] = callbayes(engine, ...
        'model'           , sprintf('%s/%s_%s.txt', modelDir, modelName, engine), ...
        'data'            , data                                      , ...
        'outputname'      , 'samples'                                 , ...
        'init'            , initArg                                   , ...
        'datafilename'    , modelName                                 , ...
        'initfilename'    , modelName                                 , ...
        'scriptfilename'  , modelName                                 , ...
        'logfilename'     , sprintf('/tmp/%s_c%d', modelName, chunk)  , ...
        'nchains'         , nChains                                   , ...
        'nburnin'         , nBurnin                                   , ...
        'nsamples'        , nTake                                     , ...
        'monitorparams'   , params                                    , ...
        'thin'            , nThin                                     , ...
        'workingdir'      , sprintf('/tmp/%s_c%d', modelName, chunk)  , ...
        'verbosity'       , 0                                         , ...
        'saveoutput'      , true                                      , ...
        'allowunderscores', 1                                         , ...
        'parallel'        , doParallel                                , ...
        'seed'            , rngSeed + 17 * chunk                      );

      if isempty(chains)
        chains = chainsChunk;
        stats = statsChunk;
        info = infoChunk;
      else
        chains = appendChains(chains, chainsChunk);
        stats = statsChunk;
        info = infoChunk;
      end

      chains = dropWrapMus(chains, nStimuli);
      [chains, nFlips] = alignMuWraps(chains, stimuli, pi, nStimuli);
      fprintf('Wrap-aligned %d chain-by-bin mu shifts\n', nFlips);

      chains = restartLostChains(chains, stimuli, pi, nStimuli, ...
        data, params, modelDir, modelName, engine, nBurninRestart, nThin, ...
        doParallel, rngSeed, chunk);

      diag = diagnosePosterior(chains, nStimuli, rHatMax, essMin, nDrawsWanted);
      printDiagnostics(diag);
      nKept = size(chains.sigma_1, 1);

      if diag.converged
        fprintf('Stop rule passed with %d draws/chain (%d joint).\n', ...
          nKept, nKept * nChains);
        break
      end
      if nKept >= nSamplesMax
        warning('memoryReproductionMCMC:noConverge', ...
          ['Reached per-chain cap (%d) without passing the stop rule. ', ...
          'Saving all %d chains anyway; inspect diagnostics.'], ...
          nSamplesMax, nChains);
        break
      end

      nTake = min(nSamplesChunk, nSamplesMax - nKept);
      nBurnin = nBurninRestart;
      setChainInitPool(lastStateInits(chains, nStimuli, nTrials, pi));
      initArg = @nextChainInit;
    end

    fprintf('%s took %f seconds!\n', upper(engine), toc);

    disp('Convergence statistics (Rhat >= 1.05):')
    grtable(chains, rHatMax)

    disp('Descriptive statistics for all chains:')
    codatable(chains);

    diagnostics = diagnosePosterior(chains, nStimuli, rHatMax, essMin, nDrawsWanted);

    fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
    if ~isfolder('storage')
      mkdir('storage');
    end
    save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
  end

  % never drop chains: wrap-align, then use all of them
  chains = dropWrapMus(chains, nStimuli);
  [chains, nFlips] = alignMuWraps(chains, stimuli, pi, nStimuli);
  if nFlips > 0
    fprintf('Wrap-aligned %d additional chain-by-bin mu shifts after load\n', nFlips);
  end
  diagnostics = diagnosePosterior(chains, nStimuli, rHatMax, essMin, nDrawsWanted);
  printDiagnostics(diagnostics);
  if ~diagnostics.converged
    warning('memoryReproductionMCMC:diagnostics', ...
      'Posterior diagnostics did not pass; summaries below still use all %d chains.', nChains);
  end

  % posterior summary for sigmas
  sigma3 = codatable(chains, 'sigma_1', @mean);
  bounds3 = prctile(chains.sigma_1(:), CI);
  fprintf('Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
  sigma6 = codatable(chains, 'sigma_2', @mean);
  bounds6 = prctile(chains.sigma_2(:), CI);
  fprintf('Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

  omega3 = get_matrix_from_coda(chains, 'omega3', @mean);
  omega6 = get_matrix_from_coda(chains, 'omega6', @mean);
  fprintf('Posterior means for omega for set size 3 are (%1.3f, %1.3f, %1.3f)\n', omega3);
  fprintf('Posterior means for omega for set size 6 are (%1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f)\n', omega6);
  for i = 1:3
    bounds = prctile(chains.(sprintf('omega3_%d', i))(:), CI);
    fprintf('Posterior mean of omega_%d for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega3(i), bounds);
  end
  for i = 1:6
    bounds = prctile(chains.(sprintf('omega6_%d', i))(:), CI);
    fprintf('Posterior mean of omega_%d for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega6(i), bounds);
  end

  % Savage-Dickey at no-swap corner, kernel density on 10k joint draws
  omega3full = [chains.omega3_1(:) chains.omega3_2(:) chains.omega3_3(:)];
  omega6full = [chains.omega6_1(:) chains.omega6_2(:) chains.omega6_3(:) ...
    chains.omega6_4(:) chains.omega6_5(:) chains.omega6_6(:)];
  omega3fullPrior = [chains.omega3prior_1(:) chains.omega3prior_2(:) chains.omega3prior_3(:)];
  omega6fullPrior = [chains.omega6prior_1(:) chains.omega6prior_2(:) chains.omega6prior_3(:) ...
    chains.omega6prior_4(:) chains.omega6prior_5(:) chains.omega6prior_6(:)];
  logBF3 = savageDickeyLogBF(omega3full, omega3fullPrior, [1 0 0]);
  logBF6 = savageDickeyLogBF(omega6full, omega6fullPrior, [1 0 0 0 0 0]);
  fprintf('Via multivariate kernel density estimation, the log BF for the null is %.0f for 3 targets and %.0f for 6 targets\n', ...
    logBF3, logBF6);
  fprintf('Joint posterior draws used for Savage-Dickey: %d\n', size(omega3full, 1));

  % inferred representation
  F = figure; clf; hold on;
  setFigure(F, [0.2 0.2 0.4 0.4], '');

  mu = nan(dm.nStimuli, 1);
  muBounds = nan(dm.nStimuli, 2);
  for idx = 1:dm.nStimuli
    vals = chains.(sprintf('mu_%d', idx))(:);
    [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CI);
  end
  muTruth = dm.stimuli;

  cla; hold on;
  set(gca, ...
    'xlim'       , [0 pi]    , ...
    'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
    'xticklabelrot', 0, ...
    'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
    'ylim'       , [0 pi]    , ...
    'ytick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
    'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
    'ticklabelinterpreter', 'latex', ...
    'box'        , 'off'     , ...
    'tickdir'    , 'out'     , ...
    'layer'      , 'top'     , ...
    'ticklength' , [0.02 0]  , ...
    'layer'      , 'top'     , ...
    'clipping'   , 'off'     , ...
    'fontsize'   , fontSize  );
  axis square;
  ylabel('Psychological', 'fontsize', fontSize);
  xlabel('Physical', 'fontsize', fontSize);
  moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
  Raxes(gca, 0.02, 0.01);

  for i = pi/4:pi/4:3*pi/4
    plot([i i], [0 pi], '-', ...
      'color', pantone.GlacierGray);
    plot([0 pi], [i i], '-', ...
      'color', pantone.GlacierGray);
  end

  for idx = 1:dm.nStimuli
    if muBounds(idx, 1) > muBounds(idx, 2)
      plot(muTruth(idx) * [1 1], [0 muBounds(idx, 2)], '-', ...
        'color', pantone.ClassicBlue, 'linewidth', 1);
      plot(muTruth(idx) * [1 1], [muBounds(idx, 1) pi], '-', ...
        'color', pantone.ClassicBlue, 'linewidth', 1);
    else
      plot(muTruth(idx) * [1 1], muBounds(idx, :), '-', ...
        'color', pantone.ClassicBlue, 'linewidth', 1);
    end
    plot(muTruth(idx), mu(idx), 'o', ...
      'markerfacecolor', pantone.ClassicBlue, ...
      'markeredgecolor', 'w', 'linewidth', 0.5, 'markersize', 4);
  end
  plot([0 pi], [0 pi], '-', ...
    'color', pantone.AuroraRed, 'linewidth', 0.5);

  if printFigures
    if ~isfolder('figures')
      mkdir('figures');
    end
    print(sprintf('figures/%s_%s.png', dataName, modelName), '-dpng');
    print(sprintf('figures/%s_%s.eps', dataName, modelName), '-depsc');
  end

end

%% local functions

function init = makeInits(stimuli, nTrials, piVal)
% Wrap-copy models declare mu[1:3*nStimuli]. Only mu[2:nStimuli] are
% stochastic; mu[1] and the +/- pi copies are deterministic — leave as NaN.
nStimuli = numel(stimuli);
mu = nan(3 * nStimuli, 1);
mu(2:nStimuli) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:nStimuli) = min(max(mu(2:nStimuli), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
omega3 = jitterSimplex([0.85 0.075 0.075]);
omega6 = jitterSimplex([0.55 0.09 0.09 0.09 0.09 0.09]);
init = struct( ...
  'mu', mu, ...
  'sigma', 0.2 + 0.7 * rand(1, 2), ...
  'omega3', omega3, ...
  'omega6', omega6, ...
  'omega3prior', jitterSimplex(ones(1, 3) / 3), ...
  'omega6prior', jitterSimplex(ones(1, 6) / 6), ...
  'xi', ones(nTrials, 2));
end

function p = jitterSimplex(p)
p = p(:)' .* exp(0.15 * randn(size(p)));
p = max(p, 1e-6);
p = p / sum(p);
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

function [chains, nFlips] = alignMuWraps(chains, stimuli, piVal, nStimuli)
nFlips = 0;
nCh = size(chains.sigma_1, 2);
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

function inits = lastStateInits(chains, nStimuli, nTrials, piVal)
nCh = size(chains.sigma_1, 2);
inits = cell(1, nCh);
for c = 1:nCh
  mu = nan(3 * nStimuli, 1);
  for j = 2:nStimuli
    fname = sprintf('mu_%d', j);
    if isfield(chains, fname)
      mu(j) = chains.(fname)(end, c);
    end
  end
  mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
  for j = 3:nStimuli
    mu(j) = min(max(mu(j), 1e-4), piVal - 1e-4);
  end
  omega3 = [chains.omega3_1(end, c) chains.omega3_2(end, c) chains.omega3_3(end, c)];
  omega3 = max(omega3, 1e-6);
  omega3 = omega3 / sum(omega3);
  omega6 = [chains.omega6_1(end, c) chains.omega6_2(end, c) chains.omega6_3(end, c) ...
    chains.omega6_4(end, c) chains.omega6_5(end, c) chains.omega6_6(end, c)];
  omega6 = max(omega6, 1e-6);
  omega6 = omega6 / sum(omega6);
  inits{c} = struct( ...
    'mu', mu, ...
    'sigma', [chains.sigma_1(end, c) chains.sigma_2(end, c)], ...
    'omega3', omega3, ...
    'omega6', omega6, ...
    'omega3prior', jitterSimplex(ones(1, 3) / 3), ...
    'omega6prior', jitterSimplex(ones(1, 6) / 6), ...
    'xi', ones(nTrials, 2));
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

function chains = restartLostChains(chains, stimuli, piVal, nStimuli, ...
  data, params, modelDir, modelName, engine, nBurnin, nThin, ...
  doParallel, rngSeed, chunk)
nCh = size(chains.sigma_1, 2);
nKept = size(chains.sigma_1, 1);
err = zeros(1, nCh);
for c = 1:nCh
  acc = 0;
  for j = 1:nStimuli
    fname = sprintf('mu_%d', j);
    if ~isfield(chains, fname)
      continue
    end
    m = circularMeanHalf(chains.(fname)(:, c), piVal);
    acc = acc + halfCircleDist(m, stimuli(j), piVal)^2;
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
nTrials = data.nTrials;
baseInits = lastStateInits(chains, nStimuli, nTrials, piVal);
template = baseInits{best};
template.sigma = max(template.sigma + 0.02 * randn(1, 2), 0.05);
template.sigma = min(template.sigma, piVal - 1e-3);
template.omega3 = jitterSimplex(template.omega3);
template.omega6 = jitterSimplex(template.omega6);
for i = 1:numel(lost)
  c = lost(i);
  [~, chainsNew, ~, ~] = callbayes(engine, ...
    'model'           , sprintf('%s/%s_%s.txt', modelDir, modelName, engine), ...
    'data'            , data                                      , ...
    'outputname'      , 'samples'                                 , ...
    'init'            , @()template                               , ...
    'datafilename'    , sprintf('%s_r%d', modelName, c)          , ...
    'initfilename'    , sprintf('%s_r%d', modelName, c)          , ...
    'scriptfilename'  , sprintf('%s_r%d', modelName, c)          , ...
    'logfilename'     , sprintf('/tmp/%s_r%d_c%d', modelName, c, chunk), ...
    'nchains'         , 1                                         , ...
    'nburnin'         , nBurnin                                   , ...
    'nsamples'        , nKept                                     , ...
    'monitorparams'   , params                                    , ...
    'thin'            , nThin                                     , ...
    'workingdir'      , sprintf('/tmp/%s_r%d_c%d', modelName, c, chunk), ...
    'verbosity'       , 0                                         , ...
    'saveoutput'      , true                                      , ...
    'allowunderscores', 1                                         , ...
    'parallel'        , doParallel                                , ...
    'seed'            , rngSeed + 100 * chunk + c                 );
  chainsNew = dropWrapMus(chainsNew, nStimuli);
  [chainsNew, ~] = alignMuWraps(chainsNew, stimuli, piVal, nStimuli);
  fn = fieldnames(chains);
  for f = 1:numel(fn)
    if isfield(chainsNew, fn{f}) && size(chainsNew.(fn{f}), 2) == 1
      chains.(fn{f})(:, c) = chainsNew.(fn{f})(:, 1);
    end
  end
end
end

function diag = diagnosePosterior(chains, nStimuli, rHatMax, essMin, nDrawsWanted)
names = monitoredNames(chains, nStimuli);
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
nKept = size(chains.sigma_1, 1);
nCh = size(chains.sigma_1, 2);
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
diag.converged = diag.maxRhat < rHatMax ...
  && diag.minEss >= essMin ...
  && diag.minEssTail >= essMin ...
  && diag.nJoint >= nDrawsWanted ...
  && nCh >= 4;
if isempty(idxR)
  diag.maxRhatName = '';
end
if isempty(idxE)
  diag.minEssName = '';
end
end

function names = monitoredNames(chains, nStimuli)
fn = fieldnames(chains);
keep = {};
for i = 1:numel(fn)
  name = fn{i};
  if strncmp(name, 'omega3prior', 11) || strncmp(name, 'omega6prior', 11)
    continue
  end
  if strncmp(name, 'mu_', 3)
    idx = sscanf(name, 'mu_%d');
    if isempty(idx) || idx < 2 || idx > nStimuli
      continue
    end
  elseif ~(strncmp(name, 'sigma_', 6) || strncmp(name, 'omega3_', 7) || strncmp(name, 'omega6_', 7))
    continue
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

function logBF = savageDickeyLogBF(posterior, priorDraws, critical)
n = size(posterior, 1);
d = size(posterior, 2);
sigma = std(posterior, 0, 1);
bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
bw(bw <= 0) = 1e-3;
fPosterior = mvksdensity(posterior, critical, 'Bandwidth', bw);
n = size(priorDraws, 1);
d = size(priorDraws, 2);
sigma = std(priorDraws, 0, 1);
bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
bw(bw <= 0) = 1e-3;
fPrior = mvksdensity(priorDraws, critical, 'Bandwidth', bw);
logBF = log(fPosterior) - log(fPrior);
end

function [muMean, bounds] = summarizeHalfCircle(vals, CI)
vals = vals(:);
phi = mod(2 * vals, 2 * pi);
phi0 = atan2(mean(sin(phi)), mean(cos(phi)));
if phi0 < 0
   phi0 = phi0 + 2 * pi;
end
rel = angle(exp(1i * (phi - phi0)));
relBounds = prctile(rel, CI);
muMean = mod(phi0 / 2, pi);
bounds = mod(muMean + relBounds / 2, pi);
end
