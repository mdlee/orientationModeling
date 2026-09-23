%% Memory reproduction with swaps (Jones-Pewsey), shared mu
%
% Doubled-angle Jones-Pewsey rewrite of memoryReproduction. Both set sizes
% share one psychological map mu; sigma, psi, and omega are set-size-
% specific. Trials are two blocks so each (kappa, psi) pair is a contiguous
% likelihood (the interleaved original hung on normalizer cache misses).
%
%   yCirc = 2 * y on [0, 2*pi)
%   yCirc3 ~ dJonesPewsey(2*mu[recalled], kappa[1], psi[1])
%   yCirc6 ~ dJonesPewsey(2*mu[recalled], kappa[2], psi[2])
%
% Requires the Jones-Pewsey JAGS module in ./jags-jonesPewsey.
% For set-size-specific mu maps see separateSetSize/withSwap.m.

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
modelName = 'memoryReproductionJP';
engine = 'jags';

dataList = {...
   'tomicBaysMemory', ...
   };

pi = 3.141592653589793;
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

         [~, ~, ssCode] = unique(dm.setSize, 'stable');
         trial3 = find(ssCode == 1);
         trial6 = find(ssCode == 2);
         nTrials3 = numel(trial3);
         nTrials6 = numel(trial6);
         nStimuli = dm.nStimuli;
         yCirc3 = mod(2 * dm.response(trial3), 2 * pi);
         yCirc6 = mod(2 * dm.response(trial6), 2 * pi);
         s3 = [dm.tIdx(trial3) dm.nIdx(trial3, 1:2)];
         s3(isnan(s3)) = 1;
         s6 = [dm.tIdx(trial6) dm.nIdx(trial6, 1:5)];
         s6(isnan(s6)) = 1;
         for i = 1:nTrials3
            t = trial3(i);
            vals = dm.nontarget(t, 1:2);
            [~, srt] = sort([0 min(abs(dm.target(t) - vals), ...
               pi - abs(dm.target(t) - vals))], 'ascend');
            s3(i, 1:3) = s3(i, srt);
         end
         for i = 1:nTrials6
            t = trial6(i);
            vals = dm.nontarget(t, 1:5);
            [~, srt] = sort([0 min(abs(dm.target(t) - vals), ...
               pi - abs(dm.target(t) - vals))], 'ascend');
            s6(i, 1:6) = s6(i, srt);
         end
         fprintf('Shared-mu swap JP: %d trials at set size 3, %d at set size 6\n', ...
            nTrials3, nTrials6);
   end

   params = {'mu', 'sigma', 'psi', 'omega3', 'omega6', 'omega3prior', 'omega6prior'};
   stimuli = dm.stimuli(:);

   data = struct(...
      's3'       , s3       , ...
      's6'       , s6       , ...
      'yCirc3'  , yCirc3   , ...
      'yCirc6'  , yCirc6   , ...
      'nStimuli', nStimuli , ...
      'nTrials3', nTrials3 , ...
      'nTrials6', nTrials6 );

   cfg = struct();
   cfg.modelFile = sprintf('%s/%s_%s.txt', modelDir, modelName, engine);
   cfg.modelName = modelName;
   cfg.data = data;
   cfg.params = params;
   cfg.nStimuli = nStimuli;
   cfg.stimuli = stimuli;
   cfg.makeInits = @()makeInits(stimuli, nTrials3, nTrials6, pi);
   cfg.lastStateInits = @(chains)lastStateInits(chains, nStimuli, nTrials3, nTrials6, pi);
   cfg.sizeField = 'sigma_1';
   cfg.modules = {'jonespewsey'};
   cfg.piVal = pi;
   cfg.restartJitter = @(t)jitterRestart(t, pi);
   cfg.extraDiagPrefixes = {'omega3', 'omega6'};

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
      warning('memoryReproductionJP:diagnostics', ...
         'Posterior diagnostics did not pass; summaries still use all %d chains.', ...
         size(chains.sigma_1, 2));
   end

   sigma3 = codatable(chains, 'sigma_1', @mean);
   bounds3 = prctile(chains.sigma_1(:), CI);
   fprintf('Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
   sigma6 = codatable(chains, 'sigma_2', @mean);
   bounds6 = prctile(chains.sigma_2(:), CI);
   fprintf('Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

   psi3 = codatable(chains, 'psi_1', @mean);
   psi3b = prctile(chains.psi_1(:), CI);
   fprintf('Posterior mean of psi for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', psi3, psi3b);
   psi6 = codatable(chains, 'psi_2', @mean);
   psi6b = prctile(chains.psi_2(:), CI);
   fprintf('Posterior mean of psi for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', psi6, psi6b);

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

function init = makeInits(stimuli, nTrials3, nTrials6, piVal)
nStimuli = numel(stimuli);
mu = nan(nStimuli, 1);
mu(2:end) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:end) = min(max(mu(2:end), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
init = struct( ...
   'mu', mu, ...
   'sigma', 0.3 + 0.4 * rand(1, 2), ...
   'psi', clipPsi(0.2 * randn(1, 2)), ...
   'omega3', jitterSimplex([0.85 0.075 0.075]), ...
   'omega6', jitterSimplex([0.55 0.09 0.09 0.09 0.09 0.09]), ...
   'omega3prior', jitterSimplex(ones(1, 3) / 3), ...
   'omega6prior', jitterSimplex(ones(1, 6) / 6), ...
   'xi3', ones(nTrials3, 1), ...
   'xi6', ones(nTrials6, 1));
end

function inits = lastStateInits(chains, nStimuli, nTrials3, nTrials6, piVal)
nCh = size(chains.sigma_1, 2);
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
   omega3 = [chains.omega3_1(end, c) chains.omega3_2(end, c) chains.omega3_3(end, c)];
   omega3 = max(omega3, 1e-6);
   omega3 = omega3 / sum(omega3);
   omega6 = [chains.omega6_1(end, c) chains.omega6_2(end, c) chains.omega6_3(end, c) ...
      chains.omega6_4(end, c) chains.omega6_5(end, c) chains.omega6_6(end, c)];
   omega6 = max(omega6, 1e-6);
   omega6 = omega6 / sum(omega6);
   inits{c} = struct( ...
      'mu', mu, ...
      'sigma', min(max([chains.sigma_1(end, c) chains.sigma_2(end, c)], 0.02), piVal - 1e-3), ...
      'psi', clipPsi([chains.psi_1(end, c) chains.psi_2(end, c)]), ...
      'omega3', omega3, ...
      'omega6', omega6, ...
      'omega3prior', jitterSimplex(ones(1, 3) / 3), ...
      'omega6prior', jitterSimplex(ones(1, 6) / 6), ...
      'xi3', ones(nTrials3, 1), ...
      'xi6', ones(nTrials6, 1));
end
end

function t = jitterRestart(t, piVal)
t.sigma = min(max(t.sigma + 0.02 * randn(size(t.sigma)), 0.05), piVal - 1e-3);
t.psi = clipPsi(t.psi + 0.05 * randn(size(t.psi)));
t.omega3 = jitterSimplex(t.omega3);
t.omega6 = jitterSimplex(t.omega6);
end

function p = jitterSimplex(p)
p = p(:)' .* exp(0.15 * randn(size(p)));
p = max(p, 1e-6);
p = p / sum(p);
end

function psi = clipPsi(psi)
psi = min(max(psi, -0.999), 0.999);
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
