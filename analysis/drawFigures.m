% figures and non-modeling analyses for tomic and bays data

clear;
close all;

printFigures = true;

analysisList = {...
   %  'threeScatter'; ...
   % 'twoRepresentations'; ...
   % 'alternativeSimilarity'; ...
   'clusterSavageDickey'; ...
   %'statisticalSummaries'; ...
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
            mu = codatable(chains, 'mu', @mean);
            muBounds = nan(ds.nStimuli, 2);
            for idx = 1:ds.nStimuli
               muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CIbounds);
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
               plot(muTruth(idx)*ones(1, 2), muBounds(idx, :),  '-', ...
                  'color', pantone.ClassicBlue, ...
                  'linewidth', 1);
               plot(muTruth(idx), mu(idx),  'o', ...
                  'markerfacecolor', pantone.ClassicBlue, ...
                  'markeredgecolor', 'w', ...
                  'linewidth', 0.5, ...
                  'markersize', 4);
            end
            plot([0 pi], [0 pi], '-', ...
               'color', pantone.AuroraRed, 'linewidth', 0.5);
         end

      case 'twoRepresentations'

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproduction_tomicBays_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays_jags'};

         fontSize = 20;

         nModels = numel(fileList);
         [nRows, nCols] = subplotArrange(nModels);

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.6 0.4], '');

         for modelIdx = 1:nModels
            fileName = fileList{modelIdx};

            fprintf('Loading pre-stored samples from file %s\n', fileName);
            load(sprintf('%s', fileName), 'chains', 'stats', 'diagnostics', 'info');

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

      case 'alternativeSimilarity'

         fileName = 'similarityComparisonAlternative_tomicBays_jags';
         fontSize = 20;
         CIbounds = [2.5 97.5];
         thresh = 2;

         fprintf('Loading pre-stored samples from file %s\n', fileName);
         load(sprintf('storage/%s', fileName), 'chains');

         % just convergent chains
         [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
         fields = fieldnames(chains);
         for i = 1:numel(fields)
            chains.(fields{i}) = chains.(fields{i})(:, keepChains);
         end

         mu = codatable(chains, 'mu', @mean);
         for idx = 1:ds.nStimuli
            muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CIbounds);
         end
         muTruth = ds.stimuli;

         F = figure; clf; hold on;
         setFigure(F, [0.2 0.2 0.6 0.5], '');

         subplot(1, 2, 1); cla; hold on;
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
         text(-1, 3.5, lower(char(64+1)), 'fontsize', fontSize, 'fontweight', 'bold');
         moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
         Raxes(gca, 0.02, 0.01);

         for i = pi/4:pi/4:3*pi/4
            plot([i i], [0 pi], '-', ...
               'color', pantone.GlacierGray);
            plot([0 pi], [i i], '-', ...
               'color', pantone.GlacierGray);
         end

         for idx = 1:ds.nStimuli
            plot( muTruth(idx)*ones(1, 2), muBounds(idx, :), '-', ...
               'color', pantone.ClassicBlue, ...
               'linewidth', 1);
            plot(muTruth(idx), mu(idx),  'o', ...
               'markerfacecolor', pantone.ClassicBlue, ...
               'markeredgecolor', 'w', ...
               'linewidth', 0.5, ...
               'markersize', 4);
         end
         plot([0 pi], [0 pi], '-', ...
            'color', pantone.AuroraRed, 'linewidth', 0.5);

         m = nan(ds.nStimuli, ds.nStimuli);
         for matchA = 1:ds.nStimuli
            for matchB = 1:ds.nStimuli
               match = find(abs(ds.aIdx - matchA) <= thresh & abs(ds.bIdx - matchB) <= thresh);
               m(matchA, matchB) = (sum(ds.response(match))+1)/(length(match)+2);
            end
         end

         subplot(1, 2, 2); cla; hold on;
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
         xlabel('Stimulus A', 'fontsize', fontSize);
         ylabel('Stimulus B', 'fontsize', fontSize);
         text(-1, 3.5, lower(char(64+2)), 'fontsize', fontSize, 'fontweight', 'bold');
         moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
         Raxes(gca, 0.01, 0.01);

         for i = 1:ds.nStimuli
            for j = 1:ds.nStimuli
               if ~isnan(m(i, j))
                  if i <= j
                     clr = m(i, j)*[1 1 1] + (1-m(i, j))*[0 0 0];
                     % clr = exp(-1/(1+(m(i, j)-1/2)/10))*pantone.Custard + (1-m(i, j))*pantone.ClassicBlue;
                  else
                     clr = [1 1 1];
                  end
                  plot((i-1)/72*pi, (j-1)/72*pi, 's', ...
                     'markerfacecolor', clr, ...
                     'markeredgecolor', 'w', ...
                     'markersize', 10);
               end
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

      case 'statisticalSummaries'

         fid = fopen('statisticalSummaries.txt', 'w');

         printFigures = false;

         % relative directory structure in github repository
         fileList = {...
            '../perceptualReproduction/storage/perceptualReproduction_tomicBays_jags'; ...
            '../memoryReproduction/storage/memoryReproduction_tomicBays_jags'; ...
            '../similarityComparison/storage/similarityComparison_tomicBays_jags'; ...
            '../clusterRepresentation/storage/clusterRepresentation_tomicBays_jags'};

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
         fprintf(fid, '\n----------\nPerceptual Reconstruction\n\n');

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
         fprintf(fid, '\n----------\nMemory Reconstruction\n\n');

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


         fclose(fid);
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