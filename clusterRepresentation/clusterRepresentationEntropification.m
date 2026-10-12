%% Cluster representation with similarity entropification
%
% Same perceptual reference, task displacements, perception, and memory
% process as clusterRepresentation. Similarity uses
% similarityComparisonEntropification (wrap-copy samples and a logistic
% in the signed distance gap, one w), centered on muS.
%
% Sampler: orientationMcmcProtocol (8 chains kept in full, data-informed
% inits, wrap alignment, split-Rhat / ESS, 10k joint draws).
%
% Storage: storage/clusterRepresentationEntropification_tomicBays_jags.mat
% Figures: figures/tomicBays_clusterRepresentationEntropification.*
% Does not read or write the original clusterRepresentation outputs.

clear; close all;
preLoad = true;
printFigures = true;

thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(fullfile(thisDir, '..', 'supportingFiles'));

modelDir = './';
modelName = 'clusterRepresentationEntropification';
engine = 'jags';

dataList = {...
   'tomicBays', ...
   };

pi = 3.1415;
load(fullfile(thisDir, '..', 'supportingFiles', 'pantoneColors.mat'), 'pantone')
fontSize = 18;
CI = [2.5 97.5];

for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBays'
         dataDir = fullfile(thisDir, '..', 'data');
         dataName = 'tomicBays';
         load(fullfile(dataDir, dataName), 'dp', 'dm', 'ds');

         nStimuli = dp.nStimuli;

         yP = dp.response;
         yM = dm.response;
         yS = ds.response;

         nPTrials = length(yP);
         sP = dp.sIdx;

         nMTrials = length(yM);
         [~, ~, setSize] = unique(dm.setSize, 'stable');
         sM = [dm.tIdx dm.nIdx];
         sM(isnan(sM)) = 1;
         [~, maxPresented] = size(sM);

         nSTrials = length(yS);
         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;
   end

   params = {'gammaPM', 'gammaPS', 'gammaPrior', ...
      'muP', 'muM', 'muS', 'delta', ...
      'sigmaP', 'sigmaM', 'sigmaS', 'w', 'omega3', 'omega6'};
   stimuli = dp.stimuli(:);

   data = struct(...
      'yP'          , yP          , ...
      'yM'          , yM          , ...
      'yS'          , yS          , ...
      'nStimuli'    , nStimuli    , ...
      'nPTrials'    , nPTrials    , ...
      'nMTrials'    , nMTrials    , ...
      'nSTrials'    , nSTrials    , ...
      'sP'          , sP          , ...
      'sM'          , sM          , ...
      'setSize'     , setSize     , ...
      'maxPresented', maxPresented, ...
      'a'           , a           , ...
      'b'           , b           , ...
      'c'           , c           , ...
      'd'           , d           );

   xAinit = clampHalfLineSample(ds.a(:), pi);
   xBinit = clampHalfLineSample(ds.b(:), pi);
   xCinit = clampHalfLineSample(ds.c(:), pi);
   xDinit = clampHalfLineSample(ds.d(:), pi);

   cfg = struct();
   cfg.modelFile = sprintf('%s/%s_%s.txt', modelDir, modelName, engine);
   cfg.modelName = modelName;
   cfg.data = data;
   cfg.params = params;
   cfg.nStimuli = nStimuli;
   cfg.stimuli = stimuli;
   cfg.makeInits = @()makeInits(stimuli, pi, nPTrials, nMTrials, nSTrials, ...
      xAinit, xBinit, xCinit, xDinit);
   cfg.lastStateInits = @(chains)lastStateInits(chains, nStimuli, pi, ...
      nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit);
   cfg.sizeField = 'sigmaS';
   cfg.piVal = pi;
   cfg.muPrefix = 'muP';
   cfg.alignCompanions = {'muM', 'muS'};
   cfg.extraDiagPrefixes = {'sigmaP', 'sigmaM', 'sigmaS', 'w', 'omega3', 'omega6', ...
      'gammaPM', 'gammaPS', 'gammaPrior', 'muM', 'muS'};
   cfg.restartJitter = @(t)jitterRestart(t, pi, nPTrials, nMTrials, nSTrials, ...
      xAinit, xBinit, xCinit, xDinit);

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
      warning('clusterRepresentationEntropification:diagnostics', ...
         'Posterior diagnostics did not pass; summaries still use all %d chains.', ...
         size(chains.sigmaS, 2));
   end

   gammaPM = codatable(chains, 'gammaPM', @mean);
   boundsPM = prctile(chains.gammaPM(:), CI);
   fprintf('Posterior mean of gammaPM is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', gammaPM, boundsPM);
   gammaPS = codatable(chains, 'gammaPS', @mean);
   boundsPS = prctile(chains.gammaPS(:), CI);
   fprintf('Posterior mean of gammaPS is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', gammaPS, boundsPS);
   gammaPrior = codatable(chains, 'gammaPrior', @mean);
   boundsPrior = prctile(chains.gammaPrior(:), CI);
   fprintf('Posterior mean of gammaPrior is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', gammaPrior, boundsPrior);

   sigmaP = codatable(chains, 'sigmaP', @mean);
   bounds = prctile(chains.sigmaP(:), CI);
   fprintf('Posterior mean of sigma perception is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigmaP, bounds);

   sigma3 = codatable(chains, 'sigmaM_1', @mean);
   bounds3 = prctile(chains.sigmaM_1(:), CI);
   fprintf('Posterior mean of sigma memory for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
   sigma6 = codatable(chains, 'sigmaM_2', @mean);
   bounds6 = prctile(chains.sigmaM_2(:), CI);
   fprintf('Posterior mean of sigma memory for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

   sigmaS = codatable(chains, 'sigmaS', @mean);
   boundsS = prctile(chains.sigmaS(:), CI);
   fprintf('Posterior mean of sigma similarity is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigmaS, boundsS);

   wMean = codatable(chains, 'w', @mean);
   wBounds = prctile(chains.w(:), CI);
   fprintf('Posterior mean of w is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', wMean, wBounds);

   gammaLo = 0; gammaHi = 3; gammaTick = 1; gammaEps = 0.052;
   gammaE = (gammaLo-gammaEps/2):gammaEps:(gammaHi + gammaEps/2);
   gammaC = gammaLo:gammaEps:gammaHi;

   countPosterior = histcounts2(chains.gammaPM(:), chains.gammaPS(:), ...
      gammaE, gammaE, ...
      'normalization', 'probability');
   count = histcounts(chains.gammaPrior(:), gammaE);
   count(1) = 2*count(1);
   count = count/sum(count);
   countPrior = count'*count;

   scale = 0.5;
   colorPrior = pantone.Custard;
   colorPosterior = pantone.ClassicBlue;

   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   set(gca, ...
      'xlim'       , [gammaLo gammaHi]    , ...
      'xtick'      , gammaLo:gammaTick:gammaHi    , ...
      'ylim'       , [gammaLo gammaHi]    , ...
      'ytick'      , gammaLo:gammaTick:gammaHi    , ...
      'box'        , 'off'     , ...
      'tickdir'    , 'out'     , ...
      'layer'      , 'top'     , ...
      'ticklength' , [0.02 0]  , ...
      'clipping'   , 'off'     , ...
      'fontsize'   , fontSize  );
   axis square;
   xlabel('$\gamma_{\mathrm{pm}}$', 'fontsize', fontSize+4, 'interp', 'latex');
   ylabel('$\gamma_{\mathrm{ps}}$', 'fontsize', fontSize+4, 'interp', 'latex');
   moveAxis(gca, [1 1 0.9 0.9], [0 0.025 0 0]);
   Raxes(gca, 0.01, 0.01);

   for i = 1:length(gammaC)
      for j = 1:length(gammaC)
         val = sqrt(countPrior(i, j))*scale;
         rectangle('position', [gammaC(i)-val/2 gammaC(j)-val/2 val val], ...
            'curvature', [1 1], ...
            'facecolor', colorPrior, ...
            'edgecolor', colorPrior);
         val = sqrt(countPosterior(i, j))*scale;
         if val > 0
            rectangle('position', [gammaC(i)-val/2 gammaC(j)-val/2 val val], ...
               'curvature', [1 1], ...
               'facecolor', colorPosterior, ...
               'edgecolor', 'w');
         end
      end
   end

   H(1) = plot(-1, -1, 'o', ...
      'markersize', 6, ...
      'markerfacecolor', colorPrior, ...
      'markeredgecolor', colorPrior);
   H(2) = plot(-1, -1, 'o', ...
      'markersize', 6, ...
      'markerfacecolor', colorPosterior, ...
      'markeredgecolor', colorPosterior);

   L = legend(H, {'prior', 'posterior'}, ...
      'fontsize', fontSize, ...
      'box', 'off', ...
      'location', 'northeast');
   set(L, 'position', get(L, 'position') + [0.075 0.05 0 0]);

   if printFigures
      if ~isfolder('figures')
         mkdir('figures');
      end
      figBase = sprintf('figures/%s_%s', dataName, modelName);
      print([figBase '.png'], '-dpng');
      print([figBase '.eps'], '-depsc');
   end

end

function init = makeInits(stimuli, piVal, nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit)
nStimuli = numel(stimuli);
muP = nan(3 * nStimuli, 1);
muP(2:nStimuli) = stimuli(2:end) + 0.05 * randn(nStimuli - 1, 1);
muP(2:nStimuli) = min(max(muP(2:nStimuli), 1e-4), piVal - 1e-4);
muP(2) = min(max(muP(2), 1e-3), piVal/2 - 1e-3);
delta = 0.1 * randn(nStimuli, 2);
delta = min(max(delta, -piVal/2 + 1e-4), piVal/2 - 1e-4);
init = struct( ...
   'muP', muP, ...
   'sigmaP', 0.1 + 0.2 * rand, ...
   'sigmaM', 0.2 + 0.7 * rand(1, 2), ...
   'sigmaS', 0.15 + 0.25 * rand, ...
   'w', 0.5 + 3.5 * rand, ...
   'gammaPM', 0.3 + 0.7 * rand, ...
   'gammaPS', 0.3 + 0.7 * rand, ...
   'gammaPrior', 0.3 + 0.7 * rand, ...
   'delta', delta, ...
   'omega3', jitterSimplex([0.85 0.075 0.075]), ...
   'omega6', jitterSimplex([0.55 0.09 0.09 0.09 0.09 0.09]), ...
   'nu', ones(nPTrials, 1), ...
   'xi', ones(nMTrials, 2), ...
   'nuS', ones(nSTrials, 1), ...
   'xA', xAinit, ...
   'xB', xBinit, ...
   'xC', xCinit, ...
   'xD', xDinit);
end

function inits = lastStateInits(chains, nStimuli, piVal, nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit)
nCh = size(chains.sigmaS, 2);
inits = cell(1, nCh);
for c = 1:nCh
   muP = nan(3 * nStimuli, 1);
   delta = zeros(nStimuli, 2);
   for j = 2:nStimuli
      fname = sprintf('muP_%d', j);
      if isfield(chains, fname)
         muP(j) = chains.(fname)(end, c);
      end
   end
   muP(2) = min(max(muP(2), 1e-3), piVal/2 - 1e-3);
   for j = 3:nStimuli
      muP(j) = min(max(muP(j), 1e-4), piVal - 1e-4);
   end
   for j = 1:nStimuli
      delta(j, 1) = chains.(sprintf('delta_%d_1', j))(end, c);
      delta(j, 2) = chains.(sprintf('delta_%d_2', j))(end, c);
   end
   delta = min(max(delta, -piVal/2 + 1e-4), piVal/2 - 1e-4);
   inits{c} = struct( ...
      'muP', muP, ...
      'sigmaP', min(max(chains.sigmaP(end, c), 1e-3), piVal - 1e-3), ...
      'sigmaM', [ ...
         min(max(chains.sigmaM_1(end, c), 1e-3), piVal - 1e-3), ...
         min(max(chains.sigmaM_2(end, c), 1e-3), piVal - 1e-3)], ...
      'sigmaS', min(max(chains.sigmaS(end, c), 1e-3), piVal - 1e-3), ...
      'w', min(max(chains.w(end, c), 1e-3), 19.9), ...
      'gammaPM', max(chains.gammaPM(end, c), 1e-3), ...
      'gammaPS', max(chains.gammaPS(end, c), 1e-3), ...
      'gammaPrior', max(chains.gammaPrior(end, c), 1e-3), ...
      'delta', delta, ...
      'omega3', jitterSimplex([ ...
         chains.omega3_1(end, c), chains.omega3_2(end, c), chains.omega3_3(end, c)]), ...
      'omega6', jitterSimplex([ ...
         chains.omega6_1(end, c), chains.omega6_2(end, c), chains.omega6_3(end, c), ...
         chains.omega6_4(end, c), chains.omega6_5(end, c), chains.omega6_6(end, c)]), ...
      'nu', ones(nPTrials, 1), ...
      'xi', ones(nMTrials, 2), ...
      'nuS', ones(nSTrials, 1), ...
      'xA', xAinit, ...
      'xB', xBinit, ...
      'xC', xCinit, ...
      'xD', xDinit);
end
end

function t = jitterRestart(t, piVal, nPTrials, nMTrials, nSTrials, xAinit, xBinit, xCinit, xDinit)
t.sigmaP = min(max(t.sigmaP + 0.02 * randn, 0.05), piVal - 1e-3);
t.sigmaM = min(max(t.sigmaM + 0.02 * randn(size(t.sigmaM)), 0.05), piVal - 1e-3);
t.sigmaS = min(max(t.sigmaS + 0.02 * randn, 0.05), piVal - 1e-3);
t.w = min(max(t.w + 0.3 * randn, 1e-3), 19.9);
t.gammaPM = max(t.gammaPM + 0.05 * randn, 1e-3);
t.gammaPS = max(t.gammaPS + 0.05 * randn, 1e-3);
t.gammaPrior = max(t.gammaPrior + 0.05 * randn, 1e-3);
t.delta = min(max(t.delta + 0.02 * randn(size(t.delta)), -piVal/2 + 1e-4), piVal/2 - 1e-4);
t.omega3 = jitterSimplex(t.omega3);
t.omega6 = jitterSimplex(t.omega6);
t.nu = ones(nPTrials, 1);
t.xi = ones(nMTrials, 2);
t.nuS = ones(nSTrials, 1);
t.xA = xAinit;
t.xB = xBinit;
t.xC = xCinit;
t.xD = xDinit;
end

function p = jitterSimplex(p)
p = p(:)' .* exp(0.15 * randn(size(p)));
p = max(p, 1e-6);
p = p / sum(p);
end

function x = clampHalfLineSample(x, piVal)
x = mod(x, piVal);
x = min(max(x, 0), piVal - 1e-10);
end
