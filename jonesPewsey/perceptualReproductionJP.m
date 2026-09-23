%% Perceptual reproduction model (Jones-Pewsey)
%
% Doubled-angle Jones-Pewsey rewrite of perceptualReproduction, with the
% shareable MCMC protocol (8 chains kept in full, data-informed inits,
% wrap alignment, split-Rhat / ESS gates, 10k joint draws).
%
%   yCirc = 2 * y on [0, 2*pi)
%   yCirc ~ dJonesPewsey(2*mu, kappa, psi)
%   psi ~ dunif(-1, 1)
%
% Requires the Jones-Pewsey JAGS module in ./jags-jonesPewsey.

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
modelName = 'perceptualReproductionJP';
engine = 'jags';

dataList = {...
   'tomicBaysPerception'; ...
   };

pi = 3.141592653589793;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysPerception'
         dataDir = fullfile(thisDir, '..', 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'dp');

         y = dp.response;
         yCirc = mod(2 * y, 2 * pi);
         s = dp.sIdx;
         nTrials = dp.nTrials;
         nStimuli = dp.nStimuli;
   end

   params = {'mu', 'sigma', 'psi'};
   stimuli = dp.stimuli(:);

   data = struct(...
      's'       , s       , ...
      'yCirc'   , yCirc   , ...
      'nStimuli', nStimuli, ...
      'nTrials' , nTrials );

   cfg = struct();
   cfg.modelFile = sprintf('%s/%s_%s.txt', modelDir, modelName, engine);
   cfg.modelName = modelName;
   cfg.data = data;
   cfg.params = params;
   cfg.nStimuli = nStimuli;
   cfg.stimuli = stimuli;
   cfg.makeInits = @()makeInits(stimuli, pi);
   cfg.lastStateInits = @(chains)lastStateInits(chains, nStimuli, pi);
   cfg.sizeField = 'sigma';
   cfg.modules = {'jonespewsey'};
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
      warning('perceptualReproductionJP:diagnostics', ...
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

   muTruth = dp.stimuli;
   mu = nan(dp.nStimuli, 1);
   muBounds = nan(dp.nStimuli, 2);
   for idx = 1:dp.nStimuli
      vals = chains.(sprintf('mu_%d', idx))(:);
      [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CI);
   end

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

   for idx = 1:dp.nStimuli
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

function init = makeInits(stimuli, piVal)
% mu[1] is deterministic (<- 0); JAGS rejects inits for non-stochastic nodes.
% Leave mu(1) as NaN so trinity/JAGS omit that element.
nStimuli = numel(stimuli);
mu = nan(nStimuli, 1);
mu(2:end) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:end) = min(max(mu(2:end), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
init = struct( ...
   'mu', mu, ...
   'sigma', 0.15 + 0.25 * rand, ...
   'psi', clipPsi(0.2 * randn));
end

function inits = lastStateInits(chains, nStimuli, piVal)
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
      'psi', clipPsi(chains.psi(end, c)));
end
end

function t = jitterRestart(t, piVal)
t.sigma = min(max(t.sigma + 0.02 * randn, 0.05), piVal - 1e-3);
t.psi = clipPsi(t.psi + 0.05 * randn);
end

function psi = clipPsi(psi)
% Keep inside (-1, 1); avoid exact endpoints where wrapped-Cauchy / cardioid
% special cases meet the Uniform(-1,1) boundary.
psi = min(max(psi, -0.999), 0.999);
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
