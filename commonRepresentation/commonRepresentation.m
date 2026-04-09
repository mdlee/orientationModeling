%% Perceptual reproduction model

clear; close all;
preLoad = true;

% graphical model script
modelDir = './';
modelName = 'commonRepresentation';
engine = 'jags';

% data sets
dataList = {...
   'tomicBays'; ...
   };


%% constants
pi = 3.1415;
load pantoneColors pantone
fontSize = 18;
CI = [2.5 97.5];

% loop over data
for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBays'
         dataDir = '../data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'dp', 'dm', 'ds');

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

          % order nontargets in s by similarity to target
         for t = 1:dm.nTrials
            vals = dm.nontarget(t, 1:(dm.setSize(t)-1));
            [dif, srt] = sort([0 min(abs(dm.target(t)-vals), pi-abs(dm.target(t)-vals))], 'ascend');
            sM(t, 1:dm.setSize(t)) = sM(t, srt);
         end

         nSTrials = length(yS);
         a = ds.aIdx;
         b = ds.bIdx;
         c = ds.cIdx;
         d = ds.dIdx;

   end

   %% sampling from graphical model
   % parameters to monitor
   params = {'mu', 'sigmaP', 'sigmaM', 'sigmaS', 'ySP', 'nu', 'xi'};

     % MCMC properties
   nChains    = 8;     % number of MCMC chains
   nBurnin    = 1e3;   % number of discarded burn-in samples
   nSamples   = 2e3;   % number of collected samples
   nThin      = 1;    % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

   % assign MATLAB variables to the observed nodes
   data = struct(...
      'yP'          , yP       , ...
      'yM'          , yM       , ...
      'yS'          , yS       , ...
      'nStimuli'    , nStimuli , ...
      'nPTrials'    , nPTrials , ...
      'nMTrials'    , nMTrials , ...
      'nSTrials'    , nSTrials , ...
      'sP'          , sP       , ...
      'sM'          , sM       , ...
      'setSize'     , setSize  , ...
      'maxPresented', maxPresented, ...
      'a'           , a        , ...
      'b'           , b        , ...
      'c'           , c        , ...
      'd'           , d        );

% censoring initial values so data have likelihood on first sample
   for t = 1:nSTrials
      if yS(t) == 0
         xAinit(t) = 0.6; xBinit(t) = 0.7;
         xCinit(t) = 0.6; xDinit(t) = 0.8;
      else
         xAinit(t) = 0.6; xBinit(t) = 0.8;
         xCinit(t) = 0.6; xDinit(t) = 0.7;
      end
   end

   % generator for initialization
   % (note intialization of mu, which encourages the
   % better log-likelihood representation)
   generator = @()struct(...
      'sigmaP', rand*pi, ...
      'xA', xAinit, ...
      'xB', xBinit, ...
      'xC', xCinit, ...
      'xD', xDinit);
   fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);

   if preLoad && isfile(sprintf('storage/%s', fileName))
      fprintf('Loading pre-stored samples for model %s on data %s\n', modelName, dataName);
      load(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info');
   else
      tic; % start clock
      [stats, chains, diagnostics, info] = callbayes(engine, ...
         'model'           , sprintf('%s/%s_%s.txt', modelDir, modelName, engine)   , ...   , ...
         'data'            , data                                      , ...
         'outputname'      , 'samples'                                 , ...
         'init'            , generator                                 , ...
         'datafilename'    , modelName                                 , ...
         'initfilename'    , modelName                                 , ...
         'scriptfilename'  , modelName                                 , ...
         'logfilename'     , sprintf('tmp/%s', modelName)              , ...
         'nchains'         , nChains                                   , ...
         'nburnin'         , nBurnin                                   , ...
         'nsamples'        , nSamples                                  , ...
         'monitorparams'   , params                                    , ...
         'thin'            , nThin                                     , ...
         'workingdir'      , sprintf('tmp/%s', modelName)              , ...
         'verbosity'       , 0                                         , ...
         'saveoutput'      , true                                      , ...
         'allowunderscores', 1                                         , ...
         'parallel'        , doParallel                                );
      fprintf('%s took %f seconds!\n', upper(engine), toc); % show timing

      % convergence of each parameter
      disp('Convergence statistics:')
      grtable(chains, 1.05)

      % basic descriptive statistics
      disp('Descriptive statistics for all chains:')
      codatable(chains);

      fprintf('Saving samples for model %s on data %s\n', modelName, dataName);
      if ~isfolder('storage')
         !mkdir storage
      end
      save(sprintf('storage/%s', fileName), 'chains', 'stats', 'diagnostics', 'info', '-v7.3');

   end

    % just convergent enough chains
   [keepChains, rHat] = findKeepChains(chains.sigmaS, 2, 1.1);
 % keepChains =  keepChains(3);
   fields = fieldnames(chains);
   for i = 1:numel(fields)
      chains.(fields{i}) = chains.(fields{i})(:, keepChains);
   end

      % posterior summary sigma
   sigmaP = codatable(chains, 'sigmaP', @mean);
   bounds = prctile(chains.sigmaP(:), CI);
   fprintf('Posterior mean of sigma perception is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigmaP, bounds);

   sigma3 = codatable(chains, 'sigmaM_1', @mean);
   bounds3 = prctile(chains.sigmaM_1(:), CI);
   fprintf('Posterior mean of sigma memory for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
   sigma6 = codatable(chains, 'sigmaM_2', @mean);
   bounds6 = prctile(chains.sigmaM_2(:), CI);
   fprintf('Posterior mean of sigma memory for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

 % inferred representation
   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   mu = nan(dp.nStimuli, 1);
   muBounds = nan(dp.nStimuli, 2);
   for idx = 1:dp.nStimuli
      vals = chains.(sprintf('mu_%d', idx))(:);
      if idx > 66 % so hard to get convergent chains here
         vals = vals(find(vals > 2));
      end
      mu(idx) = mean(vals);
      muBounds(idx, :) = prctile(vals, CI);
   end
   muTruth = dp.stimuli;

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
