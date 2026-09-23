%% Memory reproduction without swaps MCMC (wrap-copy Gaussian, automated)
%
% Same generative model as memoryReproductionNoSwap, with the shareable
% sampler: 8 chains kept in full, data-informed inits, no nu/yP monitors,
% wrap alignment, split-Rhat / ESS gates, 10k joint draws.
%
% Storage: storage/memoryReproductionNoSwapMCMC_tomicBays_jags.mat
% (does not overwrite memoryReproductionNoSwap_tomicBays_jags.mat)

clear; close all;
preLoad = true;
printFigures = true;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

modelDir = './';
modelName = 'memoryReproductionNoSwapMCMC';
engine = 'jags';

dataList = {...
   'tomicBaysMemory', ...
   };

pi = 3.1415;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysMemory'
         dataDir = fullfile(thisDir, '..', 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'dm');

         [~, ~, setSize] = unique(dm.setSize, 'stable');
         y = dm.response;
         nTrials = dm.nTrials;
         nStimuli = dm.nStimuli;
         s = dm.tIdx;
         s(isnan(s)) = 1;
   end

   params = {'mu', 'sigma'};
   stimuli = dm.stimuli(:);

   data = struct(...
      'setSize' , setSize , ...
      's'       , s       , ...
      'y'       , y       , ...
      'nStimuli', nStimuli, ...
      'nTrials' , nTrials );

   cfg = struct();
   cfg.modelFile = sprintf('%s/%s_%s.txt', modelDir, modelName, engine);
   cfg.modelName = modelName;
   cfg.data = data;
   cfg.params = params;
   cfg.nStimuli = nStimuli;
   cfg.stimuli = stimuli;
   cfg.makeInits = @()makeInits(stimuli, nTrials, pi);
   cfg.lastStateInits = @(chains)lastStateInits(chains, nStimuli, nTrials, pi);
   cfg.sizeField = 'sigma_1';
   cfg.piVal = pi;
   cfg.restartJitter = @(t)jitterRestart(t, pi);

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
      warning('memoryReproductionNoSwapMCMC:diagnostics', ...
         'Posterior diagnostics did not pass; summaries still use all %d chains.', ...
         size(chains.sigma_1, 2));
   end

   sigma3 = codatable(chains, 'sigma_1', @mean);
   bounds3 = prctile(chains.sigma_1(:), CI);
   fprintf('Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
   sigma6 = codatable(chains, 'sigma_2', @mean);
   bounds6 = prctile(chains.sigma_2(:), CI);
   fprintf('Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

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

function init = makeInits(stimuli, nTrials, piVal)
% Wrap-copy models declare mu[1:3*nStimuli]. Only mu[2:nStimuli] are
% stochastic; mu[1] and the +/- pi copies are deterministic — leave as NaN.
nStimuli = numel(stimuli);
mu = nan(3 * nStimuli, 1);
mu(2:nStimuli) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:nStimuli) = min(max(mu(2:nStimuli), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
init = struct( ...
   'mu', mu, ...
   'sigma', 0.3 + 0.4 * rand(1, 2), ...
   'nu', ones(nTrials, 1));
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
   inits{c} = struct( ...
      'mu', mu, ...
      'sigma', min(max([chains.sigma_1(end, c) chains.sigma_2(end, c)], 1e-3), piVal - 1e-3), ...
      'nu', ones(nTrials, 1));
end
end

function t = jitterRestart(t, piVal)
t.sigma = min(max(t.sigma + 0.02 * randn(size(t.sigma)), 0.05), piVal - 1e-3);
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
