%% Similarity comparison model (Jones-Pewsey)
%
% Doubled-angle Jones-Pewsey rewrite of similarityComparison, with the
% shareable MCMC protocol (8 chains kept in full, data-informed inits,
% wrap alignment, split-Rhat / ESS gates, 10k joint draws).
%
%   xAtmp..xDtmp ~ dJonesPewsey(2*mu[.], kappa, psi)
%   half-circle distance = min(|phi1-phi2|, 2*pi-|phi1-phi2|) / 2
%   y ~ dinterval(dAB - dCD, 0)
%   psi ~ dunif(-1, 1)
%
% Requires the Jones-Pewsey JAGS module in ./jags-jonesPewsey.
% Trial latents xAtmp..xDtmp must be initialized consistently with the
% censoring constraint (same helpers as the von Mises similarity driver).

clear; close all;
preLoad = true;
printFigures = true;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

jagsModuleDir = fullfile(thisDir, 'jags-jonesPewsey');
jagsModuleFile = fullfile(jagsModuleDir, 'jonespewsey.so');
if ~isfile(jagsModuleFile)
   error(['Build the JAGS Jones-Pewsey module first: cd %s && make\n' ...
      'Expected: %s'], jagsModuleDir, jagsModuleFile);
end
setenv('JAGS_LIBS', jagsModuleDir);

modelDir = './';
modelName = 'similarityComparisonJP';
engine = 'jags';

dataList = {...
   'tomicBaysSimilarity', ...
   };

pi = 3.141592653589793;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysSimilarity'
         dataDir = fullfile(thisDir, '..', 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'ds');

         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
         y = ds.response;
         nTrials = ds.nTrials;
         nStimuli = ds.nStimuli;
   end

   params = {'mu', 'sigma', 'psi'};
   stimuli = ds.stimuli(:);

   data = struct(...
      'a'        , a        , ...
      'b'        , b        , ...
      'c'        , c        , ...
      'd'        , d        , ...
      'y'        , y        , ...
      'nStimuli' , nStimuli , ...
      'nTrials'  , nTrials  );

   % Censoring inits so dinterval has positive likelihood on the first sample
   xAinit = zeros(nTrials, 1);
   xBinit = zeros(nTrials, 1);
   xCinit = zeros(nTrials, 1);
   xDinit = zeros(nTrials, 1);
   for t = 1:nTrials
      [xAinit(t), xBinit(t), xCinit(t), xDinit(t)] = censorSimilarityTrialInits( ...
         ds.a(t), ds.b(t), ds.c(t), ds.d(t), y(t), pi);
   end

   cfg = struct();
   cfg.modelFile = sprintf('%s/%s_%s.txt', modelDir, modelName, engine);
   cfg.modelName = modelName;
   cfg.data = data;
   cfg.params = params;
   cfg.nStimuli = nStimuli;
   cfg.stimuli = stimuli;
   cfg.makeInits = @()makeInits(stimuli, pi, xAinit, xBinit, xCinit, xDinit);
   cfg.lastStateInits = @(chains)lastStateInits(chains, nStimuli, pi, ...
      xAinit, xBinit, xCinit, xDinit);
   cfg.sizeField = 'sigma';
   cfg.modules = {'jonespewsey'};
   cfg.piVal = pi;
   cfg.restartJitter = @(t)jitterRestart(t, pi, xAinit, xBinit, xCinit, xDinit);

   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

   if preLoad && isfile(fullfile('storage', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(fullfile('storage', fileName), 'chains', 'stats', 'diagnostics', 'info');
      cfg.chains = chains;
      if exist('stats', 'var'), cfg.stats = stats; end
      if exist('info', 'var'), cfg.info = info; end
      [chains, stats, diagnostics, info] = orientationMcmcProtocol(cfg);
   else
      tic;
      [chains, stats, diagnostics, info] = orientationMcmcProtocol(cfg);
      fprintf('%s took %f seconds!\n', upper(engine), toc);
      fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
      if ~isfolder('storage')
         mkdir('storage');
      end
      save(fullfile('storage', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
   end

   if ~diagnostics.converged
      warning('similarityComparisonJP:diagnostics', ...
         'Posterior diagnostics did not pass; summaries still use all %d chains.', ...
         size(chains.sigma, 2));
   end

   sigma = codatable(chains, 'sigma', @mean);
   bounds = prctile(chains.sigma(:), CI);
   fprintf('Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

   psiHat = codatable(chains, 'psi', @mean);
   psiBounds = prctile(chains.psi(:), CI);
   fprintf('Posterior mean of psi is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', psiHat, psiBounds);

   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   mu = nan(ds.nStimuli, 1);
   muBounds = nan(ds.nStimuli, 2);
   for idx = 1:ds.nStimuli
      vals = chains.(sprintf('mu_%d', idx))(:);
      [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CI);
   end
   muTruth = ds.stimuli;

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
      'clipping'   , 'off'     , ...
      'fontsize'   , fontSize  );
   axis square;
   ylabel('Psychological', 'fontsize', fontSize);
   xlabel('Physical', 'fontsize', fontSize);
   moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
   Raxes(gca, 0.02, 0.01);

   for i = pi/4:pi/4:3*pi/4
      plot([i i], [0 pi], '-', 'color', pantone.GlacierGray);
      plot([0 pi], [i i], '-', 'color', pantone.GlacierGray);
   end

   for idx = 1:ds.nStimuli
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
   plot([0 pi], [0 pi], '-', 'color', pantone.AuroraRed, 'linewidth', 0.5);

   if printFigures
      if ~isfolder('figures')
         mkdir('figures');
      end
      figBase = sprintf('figures/%s_%s', dataName, modelName);
      print([figBase '.png'], '-dpng');
      print([figBase '.eps'], '-depsc');
   end

end

function init = makeInits(stimuli, piVal, xAinit, xBinit, xCinit, xDinit)
% mu[1] is deterministic (<- 0); leave as NaN so JAGS skips that element.
nStimuli = numel(stimuli);
mu = nan(nStimuli, 1);
mu(2:end) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:end) = min(max(mu(2:end), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
init = struct( ...
   'mu', mu, ...
   'sigma', 0.15 + 0.25 * rand, ...
   'psi', clipPsi(0.2 * randn), ...
   'xAtmp', xAinit, ...
   'xBtmp', xBinit, ...
   'xCtmp', xCinit, ...
   'xDtmp', xDinit);
end

function inits = lastStateInits(chains, nStimuli, piVal, xAinit, xBinit, xCinit, xDinit)
nCh = size(chains.sigma, 2);
inits = cell(1, nCh);
for c = 1:nCh
   mu = nan(nStimuli, 1);
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
   inits{c} = struct( ...
      'mu', mu, ...
      'sigma', min(max(chains.sigma(end, c), 0.02), piVal - 1e-3), ...
      'psi', clipPsi(chains.psi(end, c)), ...
      'xAtmp', xAinit, ...
      'xBtmp', xBinit, ...
      'xCtmp', xCinit, ...
      'xDtmp', xDinit);
end
end

function t = jitterRestart(t, piVal, xAinit, xBinit, xCinit, xDinit)
t.sigma = min(max(t.sigma + 0.02 * randn, 0.05), piVal - 1e-3);
t.psi = clipPsi(t.psi + 0.05 * randn);
t.xAtmp = xAinit;
t.xBtmp = xBinit;
t.xCtmp = xCinit;
t.xDtmp = xDinit;
end

function psi = clipPsi(psi)
psi = min(max(psi, -0.999), 0.999);
end

function [xA, xB, xC, xD] = censorSimilarityTrialInits(a, b, c, d, yObs, piVal)
% Minimal xAtmp inits consistent with JAGS dinterval(y).
% xAtmp are doubled-circle samples; distance matches the JP similarity model:
%   d = min(|phi1-phi2|, 2*pi-|phi1-phi2|) / 2
% y=0 => dAB < dCD; y=1 => dAB > dCD.

pi2 = 2 * piVal;
epsSep = 1e-4;
maxIter = 20;

xA = angleToDoubledSample(a, pi2);
xB = angleToDoubledSample(b, pi2);
xC = angleToDoubledSample(c, pi2);
xD = angleToDoubledSample(d, pi2);

for iter = 1:maxIter
   dAB = doubledHalfCircleDist(xA, xB, pi2);
   dCD = doubledHalfCircleDist(xC, xD, pi2);
   margin = dAB - dCD;
   ok = (yObs == 0 && margin < -epsSep) || (yObs == 1 && margin > epsSep);
   if ok
      break;
   end
   if yObs == 0
      need = max(dAB - dCD + epsSep, epsSep);
      [xA, xB] = moveDoubledPairCloser(xA, xB, need, pi2);
   else
      need = max(dCD - dAB + epsSep, epsSep);
      [xC, xD] = moveDoubledPairCloser(xC, xD, need, pi2);
   end
end

margin = doubledHalfCircleDist(xA, xB, pi2) - doubledHalfCircleDist(xC, xD, pi2);
ok = (yObs == 0 && margin < 0) || (yObs == 1 && margin > 0);
if ~ok
   error('censorSimilarityTrialInits: trial inits still inconsistent with y.');
end

xA = clampDoubled(xA, pi2);
xB = clampDoubled(xB, pi2);
xC = clampDoubled(xC, pi2);
xD = clampDoubled(xD, pi2);
end

function phi = angleToDoubledSample(angle, pi2)
phi = mod(2 * angle, pi2);
phi = min(phi, pi2 - 1e-10);
end

function phi = clampDoubled(phi, pi2)
phi = mod(phi, pi2);
phi = min(phi, pi2 - 1e-10);
end

function d = doubledHalfCircleDist(u, v, pi2)
d = min(abs(u - v), pi2 - abs(u - v)) / 2;
end

function [u, v] = moveDoubledPairCloser(u, v, halfReduction, pi2)
if halfReduction <= 0
   return;
end
circReduction = 2 * halfReduction;
diff = v - u;
if abs(diff) <= pi2 / 2
   dir = sign(diff);
   if dir == 0
      dir = 1;
   end
else
   dir = -sign(diff);
   if dir == 0
      dir = 1;
   end
end
half = circReduction / 2;
u = clampDoubled(u + dir * half, pi2);
v = clampDoubled(v - dir * half, pi2);
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
