%% Memory reproduction model with swap process

clear;
preLoad = true;
printFigures = true;

% graphical model script
modelDir = './';
modelName = 'memoryReproduction';
engine = 'jags';

% data sets
dataList = {...
  'tomicBaysMemory'; ...
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

    case 'tomicBaysMemory'
      dataDir = '../data/';
      dataName = 'tomicBays';
      load([dataDir dataName], 'dm');

      [~, ~, setSize] = unique(dm.setSize, 'stable');
      y = dm.response;
      nTrials = dm.nTrials;
      nStimuli = dm.nStimuli;
      s = [dm.tIdx dm.nIdx];
      s(isnan(s)) = 1;
      [~, maxPresented] = size(s);

      % order nontargets in s by similarity to target
      for t = 1:nTrials
        vals = dm.nontarget(t, 1:(dm.setSize(t)-1));
        [dif, srt] = sort([0 min(abs(dm.target(t)-vals), pi-abs(dm.target(t)-vals))], 'ascend');
        s(t, 1:dm.setSize(t)) = s(t, srt);
      end

  end

  %% sampling from graphical model
  % parameters to monitor
  params = {'mu', 'sigma', 'xi', 'yP', ...
    'omega3', 'omega6', 'omega3prior', 'omega6prior'};

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
    'nTrials'  , nTrials  , ...
    'maxPresented', maxPresented);

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
  keepChains = setdiff(keepChains, 8); % remove 8 via visual inspection of mu_6
  fields = fieldnames(chains);
  for i = 1:numel(fields)
    chains.(fields{i}) = chains.(fields{i})(:, keepChains);
  end



  % posterior summary for sigmas
  sigma3 = codatable(chains, 'sigma_1', @mean);
  bounds3 = prctile(chains.sigma_1(:), CI);
  fprintf('Posterior mean of sigma for set size 3 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma3, bounds3);
  sigma6 = codatable(chains, 'sigma_2', @mean);
  bounds6 = prctile(chains.sigma_2(:), CI);
  fprintf('Posterior mean of sigma for set size 6 is %1.3f, with 95%% CI (%1.3f, %1.3f)\n', sigma6, bounds6);

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


  % bayes factor to no swap model at critical points
  threshold = 0.05;
  omega3critical = [1 0 0];
  omega6critical = [1 0 0 0 0 0];

  omega3full = [chains.omega3_1(:) chains.omega3_2(:) chains.omega3_3(:)];
  [t, ~] = size(omega3full);
  diff = vecnorm(omega3full - repmat(omega3critical, t, 1), 2, 2);
  posteriorProportion3 = mean(diff < threshold);

  omega3fullPrior = [chains.omega3prior_1(:) chains.omega3prior_2(:) chains.omega3prior_3(:)];
  [t, ~] = size(omega3fullPrior);
  diff = vecnorm(omega3fullPrior - repmat(omega3critical, t, 1), 2, 2);
  priorProportion3 = mean(diff < threshold);

  omega6full = [chains.omega6_1(:) chains.omega6_2(:) chains.omega6_3(:) chains.omega6_4(:) chains.omega6_5(:) chains.omega6_6(:)];
  [t, ~] = size(omega6full);
  diff = vecnorm(omega6full - repmat(omega6critical, t, 1), 2, 2);
  posteriorProportion6 = mean(diff < threshold);

  omega6fullPrior = [chains.omega6prior_1(:) chains.omega6prior_2(:) chains.omega6prior_3(:) chains.omega6prior_4(:) chains.omega6prior_5(:) chains.omega6prior_6(:)];
  [t, ~] = size(omega6fullPrior);
  diff = vecnorm(omega6fullPrior - repmat(omega6critical, t, 1), 2, 2);
  priorProportion6 = mean(diff < threshold);

  % inferred representation
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
  end  %  % posterior predictive (to finish)
  %
  %      yP = codatable(chains, 'yP', @mean);
  %  figure
  %  plot(y, yP, 'ko');
  %
  %
  %
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

  % print
  if printFigures
    if ~isfolder('figures')
      !mkdir figures
    end
    print(sprintf('figures/%s_%s.png', dataName, modelName), '-dpng');
    print(sprintf('figures/%s_%s.eps', dataName, modelName), '-depsc');
  end

  %  % posterior predictive (to finish)
  %
  %      yP = codatable(chains, 'yP', @mean);
  %  figure
  %  plot(y, yP, 'ko');
  %
  %
  %
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
