% figures and non-modeling analyses for tomic and bays data

clear;
close all;

printFigures = true;

analysisList = {...
   % 'threeScatter'; ...
   % 'fourScatter'; ...
  % 'twoRepresentations'; ...
   %'clusterSavageDickey'; ...
   %'clusterSavageDickeyEntropification'; ...
   'statisticalSummaries'; ...
   %'descriptiveAdequacyPerceptual'; ...
   % 'descriptiveAdequacyMemory'; ...
   % 'descriptiveAdequacyMemoryNoSwap'; ...
   % 'descriptiveAdequacySimilarity'; ...
   % 'swapExamples'; ...
   };

% load data
dataDir = '../data/';
dataName = 'tomicBays';
load([dataDir dataName], 'ds', 'dp', 'dm');

% constants
load pantoneColors pantone;
pi = 3.1415;

% loops over analyses
for analysisIdx = 1:numel(analysisList)
   analysisName = analysisList{analysisIdx};

   switch analysisName
      case 'threeScatter'

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproduction_tomicBays_jags'; ...
            '../memoryReproduction/storage/memoryReproduction_tomicBays_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays_jags'};

         fontSize = 20;
         CIbounds = [2.5 97.5];
         labels = {'perceptual', 'memory', 'similarity'};

         nModels = numel(fileList);
         [nRows, nCols] = subplotArrange(nModels);

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.6 0.7], '');

         for modelIdx = 1:nModels
            fileName = fileList{modelIdx};
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

            % just keep converged chains
            switch modelIdx
               case 1 % perceptual
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
               case 2 % memory
                  [keepChains, rHat] = findKeepChains(chains.sigma_1, 2, 1.1);
                  keepChains = setdiff(keepChains, 8); % remove 8 via visual inspection of mu_6
               case 3 % similarity
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
            end
            fields = fieldnames(chains);
            for i = 1:numel(fields)
               chains.(fields{i}) = chains.(fields{i})(:, keepChains);
            end

            % posterior mean and CIs for mu
            mu = nan(ds.nStimuli, 1);
            muBounds = nan(ds.nStimuli, 2);
            for idx = 1:ds.nStimuli
               vals = chains.(sprintf('mu_%d', idx))(:);
               [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CIbounds);
            end
            muTruth = ds.stimuli;

            subplot(nRows, nCols, modelIdx); cla; hold on;
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
            if modelIdx == 1
               ylabel('Psychological', 'fontsize', fontSize);
            end
            if modelIdx == 2
               xlabel('Physical', 'fontsize', fontSize);
            end
            text(0, pi, labels{modelIdx}, ...
               'fontsize', fontSize-4, ...
               'fontweight', 'normal', ...
               'vert', 'bot', 'hor', 'lef');
            moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
            Raxes(gca, 0.02, 0.01);

            for i = pi/4:pi/4:3*pi/4
               plot([i i], [0 pi], '-', ...
                  'color', pantone.GlacierGray);
               plot([0 pi], [i i], '-', ...
                  'color', pantone.GlacierGray);
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
            plot([0 pi], [0 pi], '-', ...
               'color', pantone.AuroraRed, 'linewidth', 0.5);
         end

      case 'fourScatter'

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproduction_tomicBays_jags'; ...
            '../memoryReproduction/storage/memoryReproduction_tomicBays_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays_jags'; ...
            '../commonRepresentation/storage/commonRepresentation_tomicBays_jags'};

         fontSize = 20;
         CIbounds = [2.5 97.5];
         labels = {'perceptual', 'memory', 'similarity', 'common'};

         nModels = numel(fileList);
         [nRows, nCols] = subplotArrange(nModels);

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.6 0.7], '');

         for modelIdx = 1:nModels
            fileName = fileList{modelIdx};
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

            % just keep converged chains
            switch modelIdx
               case 1 % perceptual
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
               case 2 % memory
                  [keepChains, rHat] = findKeepChains(chains.sigma_1, 2, 1.1);
                  keepChains = setdiff(keepChains, 8); % remove 8 via visual inspection of mu_6
               case 3 % similarity
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
               case 4 % common
                  [keepChains, rHat] = findKeepChains(chains.sigmaS, 2, 1.1);
            end
            fields = fieldnames(chains);
            for i = 1:numel(fields)
               chains.(fields{i}) = chains.(fields{i})(:, keepChains);
            end

            if modelIdx ~=4
               mu = nan(dm.nStimuli, 1);
               muBounds = nan(dm.nStimuli, 2);
               for idx = 1:dm.nStimuli
                  vals = chains.(sprintf('mu_%d', idx))(:);
                  [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CIbounds);
               end
            else
               mu = nan(dp.nStimuli, 1);
               muBounds = nan(dp.nStimuli, 2);
               for idx = 1:dp.nStimuli
                  vals = chains.(sprintf('mu_%d', idx))(:);
                  if idx > 66 % so hard to get convergent chains here
                     vals = vals(vals > 2);
                  end
                  [mu(idx), muBounds(idx, :)] = summarizeHalfCircle(vals, CIbounds);
               end
            end
            muTruth = dp.stimuli;

            subplot(nRows, nCols, modelIdx); cla; hold on;
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

            text(0, pi, labels{modelIdx}, ...
               'fontsize', fontSize-4, ...
               'fontweight', 'normal', ...
               'vert', 'bot', 'hor', 'lef');
            moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
            Raxes(gca, 0.02, 0.01);

            for i = pi/4:pi/4:3*pi/4
               plot([i i], [0 pi], '-', ...
                  'color', pantone.GlacierGray);
               plot([0 pi], [i i], '-', ...
                  'color', pantone.GlacierGray);
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
            plot([0 pi], [0 pi], '-', ...
               'color', pantone.AuroraRed, 'linewidth', 0.5);
         end

         [~, T(1)] = suplabel('Physical','x');
         [~, T(2)] = suplabel('Psychological','y');
         set(T, 'fontsize', fontSize+2, 'fontweight', 'normal');
         set(T(2), 'vert', 'top');


      case 'twoRepresentations'

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproduction_tomicBays_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays_jags'; ...
            };

         fontSize = 20;

         nModels = numel(fileList);
         [nRows, nCols] = subplotArrange(nModels);

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.6 0.4], '');

         for modelIdx = 1:nModels
            fileName = fileList{modelIdx};

            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

            % just keep converged chains
            switch modelIdx
               case 1 % perceptual
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
               case 2 % similarity
                  [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
            end
            fields = fieldnames(chains);
            for i = 1:numel(fields)
               chains.(fields{i}) = chains.(fields{i})(:, keepChains);
            end

            mu = codatable(chains, 'mu', @mean);
            muTruth = ds.stimuli;

            xLim = [-1 1];
            subplot(nRows, nCols, modelIdx); cla; hold on;
            set(gca, ...
               'xlim'       , xLim    , ...
               'xtick'      , []   , ...
               'xticklabelrot', 0, ...
               'ylim'       , xLim    , ...
               'ytick'      , []  , ...
               'box'        , 'off'     , ...
               'tickdir'    , 'out'     , ...
               'layer'      , 'top'     , ...
               'ticklength' , [0.01 0]  , ...
               'layer'      , 'top'     , ...
               'clipping'   , 'off'     , ...
               'fontsize'   , fontSize  );
            moveAxis(gca, [1 1 0.9 1], [0 0 0 0]);
            axis square;
            axis off
            text(-1.15, 1.185, lower(char(64+modelIdx)), ...
               'fontsize', fontSize, 'fontweight', 'bold');

            r = [1.05 1.15];
            plot(xLim, [0 0], '-', ...
               'color', pantone.GlacierGray);
            plot(0, 0, 'o', ...
               'markerfacecolor', pantone.GlacierGray, ...
               'markeredgecolor', 'w', ...
               'markersize', 4);

            for idx = 1:ds.nStimuli
               plot([cos(muTruth(idx)) r(1)*cos(mu(idx))], [sin(muTruth(idx)) r(1)*sin(mu(idx))], '-', ...
                  'color', pantone.Titanium);
               plot(cos(muTruth(idx)), sin(muTruth(idx)), 'o', ...
                  'markerfacecolor',  'w', ...
                  'markeredgecolor', 'k', ...
                  'markersize', 5);
               plot(r*cos(mu(idx)), r*sin(mu(idx)), 'k-', ...
                  'linewidth', 2);
            end
         end


      case 'clusterSavageDickey'

         fileName = '../clusterRepresentation/storage/clusterRepresentation_tomicBays_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

         keepChains = [11 15 19]; % by visual inspection to ensure similarity follows higher-likelihood mode
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end

         % two-dimensional savage dickey analysis
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

         fontSize = 18;
         scale = 0.5;
         colorPrior = pantone.Custard;
         colorPosterior = pantone.ClassicBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.4], '');

         % axis
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

           case 'clusterSavageDickeyEntropification'

         fileName = '../clusterRepresentation/storage/clusterRepresentationEntropification_tomicBays_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

         % keepChains = [11 15 19]; % by visual inspection to ensure similarity follows higher-likelihood mode
         % fields = fieldnames(chains);
         % for i = 1:numel(fields)
         %    chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         % end

         % two-dimensional savage dickey analysis
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

         fontSize = 18;
         scale = 0.5;
         colorPrior = pantone.Custard;
         colorPosterior = pantone.ClassicBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.4], '');

         % axis
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

      case 'statisticalSummaries'

         fid = fopen('results/statisticalSummaries.txt', 'w');

         printFigures = false;

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproduction_tomicBays_jags'; ...
            '../memoryReproduction/storage/memoryReproduction_tomicBays_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays_jags'; ...
            '../clusterRepresentation/storage/clusterRepresentation_tomicBays_jags'; ...
            '../clusterRepresentation/storage/clusterRepresentationEntropification_tomicBays_jags'};

         CIbounds = [2.5 97.5];

         % PERCEPTUAL
         fileName = fileList{1};
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nPerceptual Reproduction\n\n');

         sigma = codatable(chains, 'sigma', @mean);
         bounds = prctile(chains.sigma(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

         % MEMORY
         fileName = fileList{2};
         fprintf('\n----------\nLoading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, rHat] = findKeepChains(chains.sigma_1, 2, 1.1);
         keepChains = setdiff(keepChains, 8); % remove 8 via visual inspection of mu_6
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nMemory Reproduction\n\n');

         % posterior summary for sigmas
         sigma3 = codatable(chains, 'sigma_1', @mean);
         bounds3 = prctile(chains.sigma_1(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
         sigma6 = codatable(chains, 'sigma_2', @mean);
         bounds6 = prctile(chains.sigma_2(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

         omega3 = get_matrix_from_coda(chains, 'omega3', @mean);
         omega6 = get_matrix_from_coda(chains, 'omega6', @mean);
         fprintf(fid, 'Posterior means for omega for set size 3 are (%1.3f, %1.3f, %1.3f)\n', omega3);
         fprintf(fid, 'Posterior means for omega for set size 6 are (%1.3f, %1.3f, %1.3f, %1.3f, %1.3f, %1.3f)\n', omega6);
         for i = 1:3
            bounds = prctile(chains.(sprintf('omega3_%d', i))(:), CIbounds);
            fprintf(fid, 'Posterior mean of omega_%d for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega3(i), bounds);
         end
         for i = 1:6
            bounds = prctile(chains.(sprintf('omega6_%d', i))(:), CIbounds);
            fprintf(fid, 'Posterior mean of omega_%d for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', i, omega6(i), bounds);
         end

         % bayes factor to no swap model at critical points
         omega3full = [chains.omega3_1(:) chains.omega3_2(:) chains.omega3_3(:)];
         omega6full = [chains.omega6_1(:) chains.omega6_2(:) chains.omega6_3(:) chains.omega6_4(:) chains.omega6_5(:) chains.omega6_6(:)];
         omega3fullPrior = [chains.omega3prior_1(:) chains.omega3prior_2(:) chains.omega3prior_3(:)];
         omega6fullPrior = [chains.omega6prior_1(:) chains.omega6prior_2(:) chains.omega6prior_3(:) chains.omega6prior_4(:) chains.omega6prior_5(:) chains.omega6prior_6(:)];
         omega3critical = [1 0 0];
         omega6critical = [1 0 0 0 0 0];

         % counting samples suffers from even prior having little density near critical point in the 6-dimensional case
         % threshold = 0.2;
         %
         % [t, ~] = size(omega3full);
         % diff = vecnorm(omega3full - repmat(omega3critical, t, 1), 2, 2);
         % posteriorProportion3 = mean(diff < threshold);
         %
         % [t, ~] = size(omega3fullPrior);
         % diff = vecnorm(omega3fullPrior - repmat(omega3critical, t, 1), 2, 2);
         % priorProportion3 = mean(diff < threshold);
         %
         % [t, ~] = size(omega6full);
         % diff = vecnorm(omega6full - repmat(omega6critical, t, 1), 2, 2);
         % posteriorProportion6 = mean(diff < threshold);
         %
         % [t, ~] = size(omega6fullPrior);
         % diff = vecnorm(omega6fullPrior - repmat(omega6critical, t, 1), 2, 2);
         % priorProportion6 = mean(diff < threshold);

         % based on kernel density estimation, with thanks to cursor
         n = size(omega3full, 1);
         d = size(omega3full, 2);
         sigma = std(omega3full, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = mvksdensity(omega3full, omega3critical, 'Bandwidth', bw);

         n = size(omega3fullPrior, 1);
         d = size(omega3fullPrior, 2);
         sigma = std(omega3fullPrior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = mvksdensity(omega3fullPrior, omega3critical, 'Bandwidth', bw);

         logBF3 = log(fPosterior) - log(fPrior);

         n = size(omega6full, 1);
         d = size(omega6full, 2);
         sigma = std(omega6full, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = mvksdensity(omega6full, omega6critical, 'Bandwidth', bw);

         n = size(omega6fullPrior, 1);
         d = size(omega6fullPrior, 2);
         sigma = std(omega6fullPrior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = mvksdensity(omega6fullPrior, omega6critical, 'Bandwidth', bw);

         logBF6 = log(fPosterior) - log(fPrior);

         fprintf(fid, 'Via multivariate kernel density estimation, the log BF for the null is %.0f for 3 targets and %.0f for 6 targets\n', ...
            logBF3, logBF6);

         % SIMILARITY
         fileName = fileList{3};
         fprintf('\n----------\nLoading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nSimilarity Comparison\n\n');

         sigma = codatable(chains, 'sigma', @mean);
         bounds = prctile(chains.sigma(:), CIbounds);
         fprintf(fid, 'Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

         % CLUSTER
         scale = 4;

         fileName = fileList{4};
         fprintf('\n----------\nLoading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         keepChains = [11 15 19]; % by visual inspection to ensure similarity follows higher-likelihood mode
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end
         fprintf(fid, '\n----------\nCluster Bayes Factor\n\n');

         Xposterior = [chains.gammaPS(:) chains.gammaPM(:)];
         Xprior = [chains.gammaPrior(:) chains.gammaPrior(randperm(length(chains.gammaPrior(:))))'];
         Xcrit = [0 0];
         n = size(Xposterior, 1);
         d = size(Xposterior, 2);
         sigma = std(Xposterior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = ksdensity(Xposterior, Xcrit, 'Bandwidth', bw*scale);

         n = size(Xprior, 1);
         d = size(Xprior, 2);
         sigma = std(Xprior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = ksdensity(Xprior, Xcrit, 'Bandwidth', bw*scale);

         logBF = log(fPosterior) - log(fPrior);
         fprintf(fid, 'Via kernel density estimation (with bandwidth boosted by factor of %d to avoid -inf), the log BF for the joint null is %.0f\n', ...
            scale, logBF);

         Xposterior = chains.gammaPM(:);
         Xprior = chains.gammaPrior(:);
         Xcrit = 0;
         n = size(Xposterior, 1);
         d = size(Xposterior, 2);
         sigma = std(Xposterior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = ksdensity(Xposterior, Xcrit, 'Bandwidth', bw*scale);

         n = size(Xprior, 1);
         d = size(Xprior, 2);
         sigma = std(Xprior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = ksdensity(Xprior, Xcrit, 'Bandwidth', bw*scale);

         logBF = log(fPosterior) - log(fPrior);
         fprintf(fid, 'Via kernel density estimation (with bandwidth boosted by factor of %d to avoid -inf, the log BF for the perceptual-memory null is %.0f\n', ...
            scale, logBF);

         Xposterior = chains.gammaPS(:);
         Xprior = chains.gammaPrior(:);
         Xcrit = 0;
         n = size(Xposterior, 1);
         d = size(Xposterior, 2);
         sigma = std(Xposterior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = ksdensity(Xposterior, Xcrit, 'Bandwidth', bw*scale);

         n = size(Xprior, 1);
         d = size(Xprior, 2);
         sigma = std(Xprior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = ksdensity(Xprior, Xcrit, 'Bandwidth', bw*scale);

         logBF = log(fPosterior) - log(fPrior);
         fprintf(fid, 'Via kernel density estimation (with bandwidth boosted by factor of %d to avoid -inf, the log BF for the perceptual-similarity null is %.0f\n', ...
            scale, logBF);


         % CLUSTER
         scale = 4;

         fileName = fileList{5};
         fprintf('\n----------\nLoading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');
         fprintf(fid, '\n----------\nCluster Bayes Factor (entropification\n\n');

         Xposterior = [chains.gammaPS(:) chains.gammaPM(:)];
         Xprior = [chains.gammaPrior(:) chains.gammaPrior(randperm(length(chains.gammaPrior(:))))'];
         Xcrit = [0 0];
         n = size(Xposterior, 1);
         d = size(Xposterior, 2);
         sigma = std(Xposterior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = ksdensity(Xposterior, Xcrit, 'Bandwidth', bw*scale);

         n = size(Xprior, 1);
         d = size(Xprior, 2);
         sigma = std(Xprior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = ksdensity(Xprior, Xcrit, 'Bandwidth', bw*scale);

         logBF = log(fPosterior) - log(fPrior);
         fprintf(fid, 'Via kernel density estimation (with bandwidth boosted by factor of %d to avoid -inf), the log BF for the joint null is %.0f\n', ...
            scale, logBF);

         Xposterior = chains.gammaPM(:);
         Xprior = chains.gammaPrior(:);
         Xcrit = 0;
         n = size(Xposterior, 1);
         d = size(Xposterior, 2);
         sigma = std(Xposterior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = ksdensity(Xposterior, Xcrit, 'Bandwidth', bw*scale);

         n = size(Xprior, 1);
         d = size(Xprior, 2);
         sigma = std(Xprior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = ksdensity(Xprior, Xcrit, 'Bandwidth', bw*scale);

         logBF = log(fPosterior) - log(fPrior);
         fprintf(fid, 'Via kernel density estimation (with bandwidth boosted by factor of %d to avoid -inf, the log BF for the perceptual-memory null is %.0f\n', ...
            scale, logBF);

         Xposterior = chains.gammaPS(:);
         Xprior = chains.gammaPrior(:);
         Xcrit = 0;
         n = size(Xposterior, 1);
         d = size(Xposterior, 2);
         sigma = std(Xposterior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPosterior = ksdensity(Xposterior, Xcrit, 'Bandwidth', bw*scale);

         n = size(Xprior, 1);
         d = size(Xprior, 2);
         sigma = std(Xprior, 0, 1);
         bw = sigma * (4 / (n * (d + 4))) ^ (1 / (d + 4));
         fPrior = ksdensity(Xprior, Xcrit, 'Bandwidth', bw*scale);

         logBF = log(fPosterior) - log(fPrior);
         fprintf(fid, 'Via kernel density estimation (with bandwidth boosted by factor of %d to avoid -inf, the log BF for the perceptual-similarity null is %.0f\n', ...
            scale, logBF);


         fclose(fid);


      case 'descriptiveAdequacyPerceptual'

         if exist('storage/descriptiveAdequacyPerceptualPreSave.mat', 'file')

            load('storage/descriptiveAdequacyPerceptualPreSave', 'credInterval', 'mn', 'binsC');

         else

            CI = [25 75];

            lo = 0; hi = pi; eps = 0.01;
            binsC = lo:eps:hi;
            binsE = lo-eps/2:eps:hi+eps/2;

            fileName = '../perceptualReproduction/storage/perceptualReproduction_tomicBays_jags';
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

            dp.response = mod(dp.response, pi);

            [~, ~, rIdx] = histcounts(dp.response, 'binedges', binsE);

            mn = nan(numel(binsC), 1);
            credInterval = nan(numel(binsC), 2);
            for i = 1:length(binsC)
               match = find(rIdx == i);
               nPer = numel(chains.yP_1);   % e.g. 2000*12, or use 2000*12 explicitly
               vals = zeros(numel(match) * nPer, 1);
               idx = 1;
               for j = 1:numel(match)
                  k = match(j);
                  fname = sprintf('yP_%d', k);
                  block = chains.(fname)(:);
                  vals(idx:idx+numel(block)-1) = block;
                  idx = idx + numel(block);
               end
               k = round((binsC(i) - vals) ./ pi);
               vals = vals + k .* pi;
               mn(i) = mean(vals);
               credInterval(i, :) = prctile(vals, CI);
            end

            save('storage/descriptiveAdequacyPerceptualPreSave', 'mn', 'credInterval', 'binsC');

         end

         fontSize = 18;
         scale = 0.25;
         clr = pantone.DuskBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [0 pi]    , ...
            'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
            'xticklabelrot', 0, ...
            'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ylim'       , [-0.2 pi+0.2]    , ...
            'ytick'      , [ 0 pi/4 pi/2 3*pi/4 pi ]   , ...
            'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ticklabelinterpreter', 'latex', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         % axis square;
         xlabel('Data', 'fontsize', fontSize);
         ylabel('Posterior Predictive', 'fontsize', fontSize);
         moveAxis(gca, [1 1 0.9 0.9], [0.025 0.025 0 0]);
         Raxes(gca, 0.01, 0.03);

         for i = 1:length(binsC)
            plot([binsC(i) binsC(i)], credInterval(i, :), '-', ...
               'color', pantone.Custard);
         end
         for i = 1:length(binsC)

            plot(binsC(i), mn(i), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markersize', 2, ...
               'markeredgecolor', pantone.ClassicBlue);
         end


      case 'descriptiveAdequacyMemory'

         if exist('storage/descriptiveAdequacyMemoryPreSave.mat', 'file')

            load('storage/descriptiveAdequacyMemoryPreSave', 'credInterval', 'binsC');

         else

            CI = [25 75];

            lo = 0; hi = pi; eps = 0.01;
            binsC = lo:eps:hi;
            binsE = lo-eps/2:eps:hi+eps/2;

            fileName = '../memoryReproduction/storage/memoryReproduction_tomicBays_jags';
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

            dm.response = mod(dm.response, pi);

            [~, ~, rIdx] = histcounts(dm.response, 'binedges', binsE);

            mn = nan(numel(binsC), 1);
            credInterval = nan(numel(binsC), 2);
            for i = 1:length(binsC)
               match = find(rIdx == i);
               nPer = numel(chains.yP_1);   % e.g. 2000*12, or use 2000*12 explicitly
               vals = zeros(numel(match) * nPer, 1);
               idx = 1;
               for j = 1:numel(match)
                  k = match(j);
                  fname = sprintf('yP_%d', k);
                  block = chains.(fname)(:);
                  vals(idx:idx+numel(block)-1) = block;
                  idx = idx + numel(block);
               end
               k = round((binsC(i) - vals) ./ pi);
               vals = vals + k .* pi;
               mn(i) = mean(vals);
               credInterval(i, :) = prctile(vals, CI);
            end

            save('storage/descriptiveAdequacyMemoryPreSave', 'credInterval', 'binsC');
         end

         fontSize = 18;
         scale = 0.25;
         clr = pantone.DuskBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [0 pi]    , ...
            'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
            'xticklabelrot', 0, ...
            'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ylim'       , [-0.8 pi+0.8]    , ...
            'ytick'      , [ 0 pi/4 pi/2 3*pi/4 pi ]   , ...
            'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ticklabelinterpreter', 'latex', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         xlabel('Data', 'fontsize', fontSize);
         ylabel('Posterior Predictive', 'fontsize', fontSize);
         moveAxis(gca, [1 1 0.9 0.9], [0.025 0.025 0 0]);
         Raxes(gca, 0.01, 0.03);


         for i = 1:length(binsC)
            plot([binsC(i) binsC(i)], credInterval(i, :), '-', ...
               'color', pantone.Custard);
         end
         for i = 1:length(binsC)

            plot(binsC(i), mn(i), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markersize', 2, ...
               'markeredgecolor', pantone.ClassicBlue);
         end

      case 'descriptiveAdequacyMemoryNoSwap'

         if exist('storage/descriptiveAdequacyMemoryNoSwapPreSave.mat', 'file')

            load('storage/descriptiveAdequacyMemoryNoSwapPreSave', 'credInterval', 'binsC');

         else

            CI = [25 75];

            lo = 0; hi = pi; eps = 0.01;
            binsC = lo:eps:hi;
            binsE = lo-eps/2:eps:hi+eps/2;

            fileName = '../memoryReproduction/storage/memoryReproductionNoSwap_tomicBays_jags';
            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

            dm.response = mod(dm.response, pi);

            [~, ~, rIdx] = histcounts(dm.response, 'binedges', binsE);

            mn = nan(numel(binsC), 1);
            credInterval = nan(numel(binsC), 2);
            for i = 1:length(binsC)
               match = find(rIdx == i);
               nPer = numel(chains.yP_1);   % e.g. 2000*12, or use 2000*12 explicitly
               vals = zeros(numel(match) * nPer, 1);
               idx = 1;
               for j = 1:numel(match)
                  k = match(j);
                  fname = sprintf('yP_%d', k);
                  block = chains.(fname)(:);
                  vals(idx:idx+numel(block)-1) = block;
                  idx = idx + numel(block);
               end
               k = round((binsC(i) - vals) ./ pi);
               vals = vals + k .* pi;
               mn(i) = mean(vals);
               credInterval(i, :) = prctile(vals, CI);
            end

            save('storage/descriptiveAdequacyMemoryNoSwapPreSave', 'credInterval', 'binsC');
         end

         fontSize = 18;
         scale = 0.25;
         clr = pantone.DuskBlue;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [0 pi]    , ...
            'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
            'xticklabelrot', 0, ...
            'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ylim'       , [-0.8 pi+0.8]    , ...
            'ytick'      , [ 0 pi/4 pi/2 3*pi/4 pi ]   , ...
            'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
            'ticklabelinterpreter', 'latex', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         xlabel('Data', 'fontsize', fontSize);
         ylabel('Posterior Predictive', 'fontsize', fontSize);
         moveAxis(gca, [1 1 0.9 0.9], [0.025 0.025 0 0]);
         Raxes(gca, 0.01, 0.03);


         for i = 1:length(binsC)
            plot([binsC(i) binsC(i)], credInterval(i, :), '-', ...
               'color', pantone.Custard);
         end
         for i = 1:length(binsC)

            plot(binsC(i), mn(i), 'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markersize', 2, ...
               'markeredgecolor', pantone.ClassicBlue);
         end

      case 'descriptiveAdequacySimilarity'

         nReps = 1e3;

         fileName = '../similarityComparison/storage/similarityComparison_tomicBays_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

         if exist('storage/descriptiveAdequacySimilarityPreSave.mat', 'file')

            load('storage/descriptiveAdequacySimilarityPreSave', 'pC');

         else

            nSamples = 1000*8;

            % need to sample posterior predictive because of censoring
            % base it on just posterior means of mu and sigma
            mu = zeros(ds.nStimuli, nSamples);
            sigma = chains.sigma(:);
            for i = 1:ds.nStimuli
               mu(i, :) = chains.(sprintf('mu_%d', i))(:);
            end

            muMean = codatable(chains, 'mu', @mean);
            sigmaMean = codatable(chains, 'sigma', @mean);

            pC = zeros(ds.nTrials, 1);
            for i = 1:ds.nTrials
               yP = zeros(nReps, 1);
               for rep = 1:nReps
                  xA = randn*sigma(ceil(rand*nSamples)) + mu(ds.aIdx(i), ceil(rand*nSamples));
                  xB = randn*sigma(ceil(rand*nSamples)) + mu(ds.bIdx(i), ceil(rand*nSamples));
                  xC = randn*sigma(ceil(rand*nSamples)) + mu(ds.cIdx(i), ceil(rand*nSamples));
                  xD = randn*sigma(ceil(rand*nSamples)) + mu(ds.dIdx(i), ceil(rand*nSamples));

                  dAB = min(abs(xA-xB), pi-abs(xA-xB));
                  dCD = min(abs(xC-xD), pi-abs(xC-xD));

                  yP(rep) = (dAB > dCD);
               end
               if ds.response(i) == 1
                  pC(i) = mean(yP);
               else
                  pC(i) = 1 - mean(yP);
               end
               %mean(pC(1:i))
            end

            save('storage/descriptiveAdequacySimilarityPreSave', 'pC');
         end

         fontSize = 18;
         lo = 0; hi = 1; eps = 0.02; xTick = 0.1;
         clr = pantone.ClassicBlue;

         binsC = lo:eps:hi;
         binsE = lo-eps/2:eps:hi+eps/2;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.4 0.6], '');

         set(gca, ...
            'xlim'       , [lo hi]    , ...
            'xtick'      , lo:xTick:hi   , ...
            'xticklabelrot', 0, ...
            'ycolor', 'none', ...
            'box'        , 'off'     , ...
            'tickdir'    , 'out'     , ...
            'layer'      , 'top'     , ...
            'ticklength' , [0.02 0]  , ...
            'layer'      , 'top'     , ...
            'clipping'   , 'off'     , ...
            'fontsize'   , fontSize  );
         xlabel('Posterior Predictive Agreement', 'fontsize', fontSize);
         moveAxis(gca, [1 1 1 1], [0 0.025 0 0]);
         Raxes(gca, 0.01, 0.01);

         count = histcounts(pC, 'binedges', binsE, 'norm', 'prob');
         H = bar(binsC, count);
         set(H, 'facecolor', clr, 'edgecolor', 'w', 'barwidth', 0.9);

         fprintf('\n mean agreement is %1.2f, with interquartile range %1.2f-%1.2f\n', ...
            mean(pC), prctile(pC, [2.5 97.5]));
         fprintf('\n median agreement is %d\n', ...
            round(100*median(pC)));

         fprintf('\n Forced choice agreement is %1.2f\n', mean(pC >= 0.5));

      case 'swapExamples'

         fileName = '../memoryReproduction/storage/memoryReproduction_tomicBays_jags';
         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('%s', fileName), 'chains');

         [keepChains, ~] = findKeepChains(chains.sigma_1, 2, 1.1);
         keepChains = setdiff(keepChains, 8);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end

         [T, picked] = buildSwapExampleTable(dm, chains, pi);
         if ~isfolder('figures')
            mkdir('figures');
         end
         writetable(T, 'figures/swapExampleCandidates.csv');
         writetable(picked, 'figures/swapExamplePicked.csv');

         fprintf('Visually clear candidates: %d (%d setSize=3, %d setSize=6)\n', ...
            height(T), sum(T.setSize == 3), sum(T.setSize == 6));
         fprintf('Picked %d trials (%d setSize=3, %d setSize=6)\n', ...
            height(picked), sum(picked.setSize == 3), sum(picked.setSize == 6));
         if ~isempty(picked)
            disp(picked(:, {'trial', 'setSize', 'pBestFoil', 'bestFoilXi', ...
               'dInferredFoil', 'foilMargin', 'participant'}));
            F = plotSwapExamples(picked, dm, pantone, pi);
            figure(F);
         else
            warning('swapExamples:NoTrials', ...
               'No swap-example trials passed filters; CSV tables written, figure skipped.');
         end

   end

   % print
   if printFigures
      if ~isfolder('figures')
         !mkdir figures
      end
      print(sprintf('figures/%s.png', analysisName), '-dpng');
      print(sprintf('figures/%s.eps', analysisName), '-depsc');
   end

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

function [T, picked] = buildSwapExampleTable(dm, chains, piVal)
angDist = @(a, b) min(abs(a - b), piVal - abs(a - b));
nTrials = dm.nTrials;
[~, ~, setSizeIdx] = unique(dm.setSize, 'stable');

rows = struct( ...
   'trial', {}, 'participant', {}, 'setSize', {}, ...
   'pSwap', {}, 'pTarget', {}, 'pBestFoil', {}, 'bestXi', {}, 'bestFoilXi', {}, ...
   'visualFoil', {}, 'visualMatch', {}, ...
   'resp', {}, 'target', {}, 'bestFoil', {}, 'bestFoilIdx', {}, ...
   'dTarget', {}, 'dBestFoil', {}, 'dInferredFoil', {}, ...
   'foilMargin', {}, 'score', {});

for t = 1:nTrials
   ss = setSizeIdx(t);
   xiSamples = chains.(sprintf('xi_%d_%d', t, ss))(:);
   pTarget = mean(xiSamples == 1);
   pSwap = 1 - pTarget;
   nCat = dm.setSize(t);
   counts = histcounts(xiSamples, 0.5:(nCat + 0.5));
   pCat = counts / sum(counts);
   [~, bestXi] = max(pCat);
   [pBestFoil, foilSlot] = max(pCat(2:end));
   bestFoilXi = foilSlot + 1;
   inferredFoilNum = bestFoilXi - 1;

   resp = mod(dm.response(t), piVal);
   targ = mod(dm.target(t), piVal);
   nFoils = dm.setSize(t) - 1;
   foilsRaw = mod(dm.nontarget(t, 1:nFoils), piVal);
   foilOrder = sortFoilsByTargetDistance(targ, foilsRaw, piVal);
   numberedFoils = foilsRaw(foilOrder);

   dT = angDist(resp, targ);
   dToNumbered = arrayfun(@(f) angDist(resp, numberedFoils(f)), 1:nFoils);
   [dBestFoil, visualFoil] = min(dToNumbered);
   foilMargin = dT - dBestFoil;
   dInferredFoil = dToNumbered(inferredFoilNum);
   visualMatch = visualFoil == inferredFoilNum;

   if dm.setSize(t) == 3
      minMargin = 0.8;
      maxDInferred = 0.35;
      minPFoil = 0.33;
   elseif bestFoilXi ~= 6
      minMargin = 0.85;
      maxDInferred = 0.05;
      minPFoil = 0.20;
   else
      minMargin = 1.0;
      maxDInferred = 0.30;
      minPFoil = 0.35;
   end

   if foilMargin < minMargin || dInferredFoil > maxDInferred || pBestFoil < minPFoil || ~visualMatch
      continue;
   end

   score = pBestFoil * foilMargin / (dInferredFoil + 0.05);
   rows(end+1) = struct( ...
      'trial', t, ...
      'participant', dm.participant(t), ...
      'setSize', dm.setSize(t), ...
      'pSwap', pSwap, ...
      'pTarget', pTarget, ...
      'pBestFoil', pBestFoil, ...
      'bestXi', bestXi, ...
      'bestFoilXi', bestFoilXi, ...
      'visualFoil', visualFoil, ...
      'visualMatch', visualMatch, ...
      'resp', resp, ...
      'target', targ, ...
      'bestFoil', numberedFoils(visualFoil), ...
      'bestFoilIdx', foilOrder(visualFoil), ...
      'dTarget', dT, ...
      'dBestFoil', dBestFoil, ...
      'dInferredFoil', dInferredFoil, ...
      'foilMargin', foilMargin, ...
      'score', score); %#ok<AGROW>
end

T = struct2table(rows);
T = sortrows(T, 'score', 'descend');

nPick = 6;
T3 = T(T.setSize == 3, :);
T6 = T(T.setSize == 6, :);
pick3 = pickDiverseSwapExamples(T3, nPick);
pick6 = pickSixItemSwapExamples(T6, nPick);
picked = [pick3; pick6];
picked = sortrows(picked, {'setSize', 'score'}, {'ascend', 'descend'});
end

function sub = pickDiverseSwapExamples(T, nPick)
if height(T) < nPick
   sub = T;
   return;
end
sub = T(1, :);
usedP = sub.participant(1);
rem = T(2:end, :);
while height(sub) < nPick && ~isempty(rem)
   best = 1;
   bestScore = -inf;
   for i = 1:height(rem)
      p = rem.participant(i);
      bonus = 0.15 * ~ismember(p, usedP);
      s = rem.score(i) * (1 + bonus);
      if s > bestScore
         bestScore = s;
         best = i;
      end
   end
   sub = [sub; rem(best, :)]; %#ok<AGROW>
   usedP(end+1) = rem.participant(best);
   rem(best, :) = [];
end
end

function pick6 = pickSixItemSwapExamples(T6, nPick)
Tother = T6(T6.bestFoilXi ~= 6, :);
Tother = sortrows(Tother, {'dInferredFoil', 'foilMargin'}, {'ascend', 'descend'});
Tfive = sortrows(T6(T6.bestFoilXi == 6, :), 'score', 'descend');

if isempty(Tother)
   pick6 = pickDiverseSwapExamples(Tfive, nPick);
   return;
end

pickOther = pickDiverseSwapExamples(Tother, 1);
usedP = pickOther.participant;
TfiveRem = Tfive(~ismember(Tfive.participant, usedP), :);
pickFive = pickDiverseSwapExamples(TfiveRem, nPick - 1);
pick6 = [pickOther; pickFive];
end

function fig = plotSwapExamples(picked, dm, pantone, piVal)
n = height(picked);
if n < 1
   error('plotSwapExamples:EmptyTable', 'picked must contain at least one trial.');
end
nCols = 4;
nRows = ceil(n / nCols);
fig = figure('Color', 'w', 'Position', [50 50 820 560], 'Units', 'pixels');
tlo = tiledlayout(nRows, nCols, 'TileSpacing', 'none', 'Padding', 'none');

xLim = [-1.12 1.12];
yLim = [-0.02 1.06];
rayLen = 0.94;
labelFs = 11;
legendFs = 12;
legendHandles = gobjects(3, 1);

for k = 1:n
   tr = picked.trial(k);
   ax = nexttile(tlo);
   hold(ax, 'on');

   targ = mod(dm.target(tr), piVal);
   nFoils = dm.setSize(tr) - 1;
   foilsRaw = mod(dm.nontarget(tr, 1:nFoils), piVal);
   resp = mod(dm.response(tr), piVal);
   foilOrder = sortFoilsByTargetDistance(targ, foilsRaw, piVal);
   inferredFoil = picked.bestFoilXi(k) - 1;

   setupSemicircleAxes(ax, xLim, yLim, pantone);
   drawSemicircleRay(ax, targ, pantone.ClassicBlue, '-', 2.8, rayLen, 0.45);
   drawSemicircleRay(ax, resp, pantone.AuroraRed, '-', 2.8, rayLen, 0.45);
   hFoil = gobjects(1);
   for fIdx = 1:nFoils
      ang = foilsRaw(foilOrder(fIdx));
      isInferred = fIdx == inferredFoil;
      h = drawSemicircleRay(ax, ang, [0 0 0], '--', 1.4, rayLen * 0.92);
      if fIdx == 1
         hFoil = h;
      end
      labelFoilNumber(ax, ang, fIdx, rayLen, isInferred, labelFs);
   end

   if k == 1
      legendHandles(1) = plot(ax, NaN, NaN, '-', 'Color', pantone.ClassicBlue, 'LineWidth', 2.8);
      legendHandles(2) = hFoil;
      legendHandles(3) = plot(ax, NaN, NaN, '-', 'Color', pantone.AuroraRed, 'LineWidth', 2.8);
   end
end

for k = (n + 1):(nRows * nCols)
   nexttile(tlo);
   axis off;
end

tlo.InnerPosition(1) = 0.02;
tlo.InnerPosition(3) = 0.96;
tlo.InnerPosition(2) = 0.04;
tlo.InnerPosition(4) = 0.90;

lgd = legend(legendHandles, {'Target', 'Foil', 'Response'}, ...
   'Orientation', 'horizontal', 'Box', 'off', 'FontSize', legendFs);
lgd.Units = 'normalized';
lgd.Position = [0.24 0.945 0.52 0.045];
end

function foilOrder = sortFoilsByTargetDistance(targ, foils, piVal)
dists = arrayfun(@(f) min(abs(targ - f), piVal - abs(targ - f)), foils);
[~, foilOrder] = sort(dists, 'ascend');
end

function labelFoilNumber(ax, angle, foilNum, rayLen, circleInferred, fontSize)
labelR = rayLen * 1.07;
x = labelR * cos(angle);
y = labelR * sin(angle);
text(ax, x, y, sprintf('%d', foilNum), ...
   'HorizontalAlignment', 'center', ...
   'VerticalAlignment', 'middle', ...
   'FontSize', fontSize, ...
   'Color', [0 0 0], ...
   'Clipping', 'off');
if circleInferred
   r = 0.050 + 0.0035 * fontSize;
   pos = [x - r, y - r, 2 * r, 2 * r];
   rectangle(ax, 'Position', pos, ...
      'Curvature', [1 1], ...
      'EdgeColor', [0 0 0], ...
      'LineWidth', 1.1, ...
      'FaceColor', 'none');
end
end

function setupSemicircleAxes(ax, xLim, yLim, pantone)
set(ax, ...
   'xlim', xLim, ...
   'ylim', yLim, ...
   'xtick', [], ...
   'ytick', [], ...
   'box', 'off', ...
   'clipping', 'off', ...
   'layer', 'top');
axis(ax, 'equal');
axis(ax, 'off');

th = linspace(0, pi, 120);
plot(ax, cos(th), sin(th), '-', 'Color', pantone.GlacierGray, 'LineWidth', 0.9);
plot(ax, xLim, [0 0], '-', 'Color', pantone.GlacierGray, 'LineWidth', 0.9);
plot(ax, 0, 0, 'o', ...
   'markerfacecolor', pantone.GlacierGray, ...
   'markeredgecolor', 'w', ...
   'markersize', 2.5);
end

function h = drawSemicircleRay(ax, angle, col, ls, lw, len, alpha)
if nargin < 7
   alpha = 1;
end
angle = mod(angle, pi);
if numel(col) == 3
   col = [col, alpha];
end
h = plot(ax, [0 len * cos(angle)], [0 len * sin(angle)], ...
   'Color', col, 'LineStyle', ls, 'LineWidth', lw);
end
