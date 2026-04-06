%% Memory reproduction model with swap process

clear;
preLoad = true;
printFigures = true;

% graphical model script
modelDir = './';
modelName = 'memoryReproductionNoSwap';
engine = 'jags';
CI = [2.5 97.5];

% data sets
dataList = {...
   'tomicBaysMemory'; ...
   };

%% general constants
pi = 3.1415;

% loop over data
for dataIdx = 1:numel(dataList)
   dataName = dataList{dataIdx};
   switch dataName

      case 'tomicBaysMemory'
         dataDir = '../data/';
         dataName = 'tomicBays';
         load([dataDir dataName], 'dm');
 
         [~, ~, setSize] = unique(dm.setSize, 'stable');
         y = dm.response;
         nTrials = dm.nTrials;
         nStimuli = dm.nStimuli;
         s = dm.tIdx;
         s(isnan(s)) = 1;
         [~, maxPresented] = size(s);

   end

   %% sampling from graphical model
   % parameters to monitor
   params = {'mu', 'sigma', 'nu', 'yP', 'omega3', 'omega6'};

   % MCMC properties
   % MCMC properties
   nChains    = 12;     % number of MCMC chains
   nBurnin    = 2e3;   % number of discarded burn-in samples
   nSamples   = 2e3;   % number of collected samples
   nThin      = 5;    % number of samples between those collected
   doParallel = 1;     % whether MATLAB parallel toolbox parallizes chains

   % assign MATLAB variables to the observed nodes
 data = struct(...
       'setSize'  , setSize  , ...
     's'        , s        , ...
      'y'        , y        , ...
      'nStimuli' , nStimuli , ...
      'nTrials'  , nTrials  );

   % generator for initialization
   generator = @()struct('sigma', rand(1, 2)*pi);

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
   [keepChains, rHat] = findKeepChains(chains.sigma_1, 2, 1.1);
   fields = fieldnames(chains);
   for i = 1:numel(fields)
      chains.(fields{i}) = chains.(fields{i})(:, keepChains);
   end

      yP = codatable(chains, 'yP', @mean);
   figure
   plot(y, yP, 'ko');


   % posterior summary for sigmas
   sigma3 = codatable(chains, 'sigma_1', @mean);
   bounds3 = prctile(chains.sigma_1(:), CI);
   fprintf('Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
   sigma6 = codatable(chains, 'sigma_2', @mean);
   bounds6 = prctile(chains.sigma_2(:), CI);
   fprintf('Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);
  
% figures

% inferred representation

   fontSize = 18;
   load pantoneColors pantone
   CIbounds = [2.5 97.5];

   F = figure; clf; hold on;
   setFigure(F, [0.2 0.2 0.4 0.4], '');

   mu = codatable(chains, 'mu', @mean);
   muBounds = nan(dm.nStimuli, 2);
   for idx = 1:dm.nStimuli
      muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CI);
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

   for idx = 1:dm.nStimuli
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

  %  % posterior predictive (to finish)
  %  lo = 0; hi = pi; step = 1/72*pi; overflow = 20;
  %  binsC = lo-overflow*step:step:hi+overflow*step;
  %  binsE = low-overflow*step-step/2:step:hi+overflow*step+step/2;
  % 
  %  F = figure; clf; hold on;
  %  setFigure(F, [0.2 0.2 0.4 0.4], '');
  % 
  %  mu = codatable(chains, 'mu', @mean);
  %  muBounds = nan(dm.nStimuli, 2);
  %  for idx = 1:dm.nStimuli
  %     muBounds(idx, :) = prctile(chains.(sprintf('mu_%d', idx))(:), CIbounds);
  %  end
  %  muTruth = dm.stimuli;
  % 
  % cla; hold on;
  %  set(gca, ...
  %     'xlim'       , [0 pi]    , ...
  %     'xtick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
  %     'xticklabelrot', 0, ...
  %     'xticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
  %     'ylim'       , [0 pi]    , ...
  %     'ytick'      , [0 pi/4 pi/2 3*pi/4 pi]   , ...
  %     'yticklabel' , {'$0$', '$\frac{\pi}{4}$', '$\frac{\pi}{2}$', '$\frac{3\pi}{4}$', '$\pi$'}, ...
  %     'ticklabelinterpreter', 'latex', ...
  %     'box'        , 'off'     , ...
  %     'tickdir'    , 'out'     , ...
  %     'layer'      , 'top'     , ...
  %     'ticklength' , [0.02 0]  , ...
  %     'layer'      , 'top'     , ...
  %     'clipping'   , 'off'     , ...
  %     'fontsize'   , fontSize  );
  %  axis square;
  %  ylabel('Psychological', 'fontsize', fontSize);
  %  xlabel('Physical', 'fontsize', fontSize);
  %  moveAxis(gca, [1 1 0.95 0.95], [0 0.025 0 0]);
  %  Raxes(gca, 0.02, 0.01);
  % 
  %  for i = pi/4:pi/4:3*pi/4
  %     plot([i i], [0 pi], '-', ...
  %        'color', pantone.GlacierGray);
  %     plot([0 pi], [i i], '-', ...
  %        'color', pantone.GlacierGray);
  %  end
  % 
  %  for idx = 1:dm.nStimuli
  %     plot(muTruth(idx)*ones(1, 2), muBounds(idx, :),  '-', ...
  %        'color', pantone.ClassicBlue, ...
  %        'linewidth', 1);
  %     plot(muTruth(idx), mu(idx),  'o', ...
  %        'markerfacecolor', pantone.ClassicBlue, ...
  %        'markeredgecolor', 'w', ...
  %        'linewidth', 0.5, ...
  %        'markersize', 4);
  % 
  %  end
  %  plot([0 pi], [0 pi], '-', ...
  %     'color', pantone.AuroraRed, 'linewidth', 0.5);

end
