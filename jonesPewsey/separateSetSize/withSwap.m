%% With-swap memory reproduction, one Jones-Pewsey fit per set size
%
% Each set size (3 and 6) is a separate analysis with its own mu, sigma, psi,
% and omega. This is not the paper-style shared-mu model.
%
%   yCirc = 2 * y on [0, 2*pi)
%   recalled item xi mixes over similarity-ranked nontargets
%   yCirc ~ dJonesPewsey(2*mu[recalled], kappa, psi)
%   psi ~ dunif(-1, 1)
%
% Requires the Jones-Pewsey JAGS module in ../jags-jonesPewsey.
% Shared-mu (both set sizes in one model) lives in the parent folder:
%   ../memoryReproductionJP.m

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
   'tomicBaysMemory', ...
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
      modelName = sprintf('withSwap_setSize%d', ss);
      trialIdx = find(ssCode == find(setSizes == ss));
      nTrials = numel(trialIdx);
      if nTrials < 1
         warning('memoryReproductionJP:noTrials', ...
            'No trials with set size %d; skipping.', ss);
         continue
      end

      yCirc = mod(2 * dm.response(trialIdx), 2 * pi);
      s = [dm.tIdx(trialIdx) dm.nIdx(trialIdx, 1:(ss - 1))];
      s(isnan(s)) = 1;
      for i = 1:nTrials
         t = trialIdx(i);
         vals = dm.nontarget(t, 1:(ss - 1));
         [~, srt] = sort([0 min(abs(dm.target(t) - vals), ...
            pi - abs(dm.target(t) - vals))], 'ascend');
         s(i, 1:ss) = s(i, srt);
      end

      fprintf('Set size %d: %d trials\n', ss, nTrials);

      params = {'mu', 'sigma', 'psi', 'omega', 'omegaPrior'};
      data = struct(...
         's'              , s                    , ...
         'yCirc'          , yCirc                , ...
         'nStimuli'       , nStimuli              , ...
         'nTrials'        , nTrials               , ...
         'omegaPriorProbs', ones(1, ss)           );

      cfg = struct();
      cfg.modelFile = fullfile(thisDir, 'withSwap_jags.txt');
      cfg.modelName = modelName;
      cfg.data = data;
      cfg.params = params;
      cfg.nStimuli = nStimuli;
      cfg.stimuli = stimuli;
      cfg.makeInits = @()makeInits(stimuli, nTrials, ss, pi);
      cfg.lastStateInits = @(chains)lastStateInits(chains, nStimuli, ss, nTrials, pi);
      cfg.sizeField = 'sigma';
      cfg.modules = {'jonespewsey'};
      cfg.piVal = pi;
      cfg.restartJitter = @(t)jitterRestart(t, ss, pi);
      cfg.extraDiagPrefixes = {'omega'};

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
         warning('memoryReproductionJP:diagnostics', ...
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

      omegaHat = get_matrix_from_coda(chains, 'omega', @mean);
      fprintf('Set size %d: posterior means for omega are (', ss);
      fprintf('%1.3f', omegaHat(1));
      fprintf(', %1.3f', omegaHat(2:end));
      fprintf(')\n');
      for i = 1:ss
         bounds = prctile(chains.(sprintf('omega_%d', i))(:), CI);
         fprintf('Set size %d: posterior mean of omega_%d is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', ...
            ss, i, omegaHat(i), bounds);
      end

      omegaFull = zeros(numel(chains.omega_1(:)), ss);
      omegaFullPrior = zeros(numel(chains.omegaPrior_1(:)), ss);
      for i = 1:ss
         omegaFull(:, i) = chains.(sprintf('omega_%d', i))(:);
         omegaFullPrior(:, i) = chains.(sprintf('omegaPrior_%d', i))(:);
      end
      critical = [1 zeros(1, ss - 1)];
      logBF = savageDickeyLogBF(omegaFull, omegaFullPrior, critical);
      fprintf('Set size %d: log BF for the null (via multivariate KDE) is %.0f\n', ss, logBF);
      fprintf('Set size %d: joint posterior draws used for Savage-Dickey: %d\n', ...
         ss, size(omegaFull, 1));

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

function init = makeInits(stimuli, nTrials, nPresented, piVal)
% mu[1] is deterministic (<- 0); leave as NaN so JAGS skips that element.
nStimuli = numel(stimuli);
mu = nan(nStimuli, 1);
mu(2:end) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
mu(2:end) = min(max(mu(2:end), 1e-4), piVal - 1e-4);
mu(2) = min(max(mu(2), 1e-3), piVal/2 - 1e-3);
if nPresented == 3
   omega0 = [0.85 0.075 0.075];
else
   omega0 = [0.55 0.09 0.09 0.09 0.09 0.09];
end
init = struct( ...
   'mu', mu, ...
   'sigma', 0.3 + 0.4 * rand, ...
   'psi', clipPsi(0.2 * randn), ...
   'omega', jitterSimplex(omega0), ...
   'omegaPrior', jitterSimplex(ones(1, nPresented) / nPresented), ...
   'xi', ones(nTrials, 1));
end

function inits = lastStateInits(chains, nStimuli, nPresented, nTrials, piVal)
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
   omega = zeros(1, nPresented);
   for i = 1:nPresented
      omega(i) = chains.(sprintf('omega_%d', i))(end, c);
   end
   omega = max(omega, 1e-6);
   omega = omega / sum(omega);
   inits{c} = struct( ...
      'mu', mu, ...
      'sigma', min(max(chains.sigma(end, c), 0.02), piVal - 1e-3), ...
      'psi', clipPsi(chains.psi(end, c)), ...
      'omega', omega, ...
      'omegaPrior', jitterSimplex(ones(1, nPresented) / nPresented), ...
      'xi', ones(nTrials, 1));
end
end

function t = jitterRestart(t, nPresented, piVal)
t.sigma = min(max(t.sigma + 0.02 * randn, 0.05), piVal - 1e-3);
t.psi = clipPsi(t.psi + 0.05 * randn);
t.omega = jitterSimplex(t.omega);
t.omegaPrior = jitterSimplex(ones(1, nPresented) / nPresented);
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
