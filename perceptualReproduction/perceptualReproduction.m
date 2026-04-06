%% Perceptual reproduction model

clear; close all;
preLoad = true;
printFigures = true;

% graphical model script
modelDir = './';
modelName = 'perceptualReproduction';
engine = 'jags';

% data sets
dataList = {...
   'tomicBaysPerception'; ...
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

      case 'tomicBaysPerception'
         dataDir = '../data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'dp');

         y = dp.response;
         s = dp.sIdx;
         nTrials = dp.nTrials;
         nStimuli = dp.nStimuli;
   end

   %% sampling from graphical model
   % parameters to monitor
   params = {'mu', 'sigma', 'nu', 'yP'};

   % MCMC properties
   nChains    = 12;     % number of MCMC chains
   nBurnin    = 2e3;   % number of discarded burn-in samples
   nSamples   = 2e3;   % number of collected samples
   nThin      = 5;    % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains
 
   % assign MATLAB variables to the observed nodes
   data = struct(...
      's'       , s       , ...
      'y'       , y       , ...
      'nStimuli', nStimuli, ...
      'nTrials' , nTrials );

   % generator for initialization
   generator = @()struct('sigma', rand*pi);

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
         'logfilename'     , sprintf('/tmp/%s', modelName)              , ...
         'nchains'         , nChains                                   , ...
         'nburnin'         , nBurnin                                   , ...
         'nsamples'        , nSamples                                  , ...
         'monitorparams'   , params                                    , ...
         'thin'            , nThin                                     , ...
         'workingdir'      , sprintf('/tmp/%s', modelName)              , ...
         'verbosity'       , 0                                         , ...
         'saveoutput'      , true                                      , ...
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
   [keepChains, rHat] = findKeepChains(chains.sigma, 2, 1.1);
   fields = fieldnames(chains);
   for i = 1:numel(fields)
      chains.(fields{i}) = chains.(fields{i})(:, keepChains);
   end

   % posterior summary sigma
   sigma = codatable(chains, 'sigma', @mean);
   bounds = prctile(chains.sigma(:), CI);
   fprintf('Posterior mean of sigma is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma, bounds);

   % inferred representation
   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   mu = codatable(chains, 'mu', @mean);
   muBounds = nan(dp.nStimuli, 2);
   for idx = 1:dp.nStimuli
      muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CI);
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

    % print
   if printFigures 
      if ~isfolder('figures')
         !mkdir figures
      end
      print(sprintf('figures/%s_%s.png', dataName, modelName), '-dpng');
      print(sprintf('figures/%s_%s.eps', dataName, modelName), '-depsc');
   end

end
