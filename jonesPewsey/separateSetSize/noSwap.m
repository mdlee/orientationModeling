%% No-swap memory reproduction, one Jones-Pewsey fit per set size
%
% Each set size (3 and 6) is a separate analysis with its own mu, sigma, and
% psi. This is not the paper-style shared-mu model.
%
%   yCirc = 2 * y on [0, 2*pi)
%   yCirc ~ dJonesPewsey(2*mu, kappa, psi)
%   psi ~ dunif(-1, 1)
%
% Requires the Jones-Pewsey JAGS module in ../jags-jonesPewsey.
% Shared-mu (both set sizes in one model) lives in the parent folder:
%   ../memoryReproductionNoSwapJP.m

clear; close all;
preLoad = true;
printFigures = true;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
repoRoot = fullfile(thisDir, '..', '..');
addpath(fullfile(repoRoot, 'supportingFiles'));

jagsModuleDir = fullfile(thisDir, '..', 'jags-jonesPewsey');
jagsModuleFile = fullfile(jagsModuleDir, 'jonespewsey.so');
if ~isfile(jagsModuleFile)
   error(['Build the JAGS Jones-Pewsey module first: cd %s && make\n' ...
      'Expected: %s'], jagsModuleDir, jagsModuleFile);
end
setenv('JAGS_LIBS', jagsModuleDir);

engine = 'jags';

dataList = {...
   'tomicBaysMemory'; ...
   };

pi = 3.141592653589793;
load(fullfile(repoRoot, 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];
setSizes = [3 6];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysMemory'
         dataDir = fullfile(repoRoot, 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'dm');
   end

   stimuli = dm.stimuli(:);
   nStimuli = dm.nStimuli;
   [~, ~, ssCode] = unique(dm.setSize, 'stable');

   for ss = setSizes
      modelName = sprintf('noSwap_setSize%d', ss);
      trialIdx = find(ssCode == find(setSizes == ss));
      nTrials = numel(trialIdx);
      if nTrials < 1
         warning('memoryReproductionNoSwapJP:noTrials', ...
            'No trials with set size %d; skipping.', ss);
         continue
      end

      yCirc = mod(2 * dm.response(trialIdx), 2 * pi);
      s = dm.tIdx(trialIdx);
      s(isnan(s)) = 1;

      fprintf('Set size %d: %d trials\n', ss, nTrials);

      params = {'mu', 'sigma', 'psi'};
      data = struct(...
         's'       , s       , ...
         'yCirc'   , yCirc   , ...
         'nStimuli', nStimuli, ...
         'nTrials' , nTrials );

      cfg = struct();
      cfg.modelFile = fullfile(thisDir, 'noSwap_jags.txt');
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
         fprintf('%s (set size %d) took %f seconds!\n', upper(engine), ss, toc);
         fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
         if ~isfolder('storage')
            mkdir('storage');
         end
         save(fullfile('storage', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
      end

      if ~diagnostics.converged
         warning('memoryReproductionNoSwapJP:diagnostics', ...
            'Set size %d: posterior diagnostics did not pass; summaries still use all %d chains.', ...
            ss, size(chains.sigma, 2));
      end

      sigmaHat = codatable(chains, 'sigma', @mean);
      sigmaB = prctile(chains.sigma(:), CI);
      fprintf('Set size %d: posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', ...
         ss, sigmaHat, sigmaB);

      psiHat = codatable(chains, 'psi', @mean);
      psiB = prctile(chains.psi(:), CI);
      fprintf('Set size %d: posterior mean of psi is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', ...
         ss, psiHat, psiB);

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
end

function init = makeInits(stimuli, piVal)
% mu[1] is deterministic (<- 0); JAGS rejects inits for non-stochastic nodes.
nStimuli = numel(stimuli);
mu = nan(nStimuli, 1);
mu(2:end) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:end) = min(max(mu(2:end), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
init = struct( ...
   'mu', mu, ...
   'sigma', 0.3 + 0.4 * rand, ...
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
