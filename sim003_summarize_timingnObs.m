%% Summarize simulation 003: timming as a function of number of obserations
% Set paths and initialize
rootDir         = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation003_CompareTimes_nObs';
settings        = load(fullfile(rootDir, 'Settings.mat'), 'settings');
settings        = settings.settings;
header          = {'NumObservations', 'FEMA_FE', 'FEMA_SE', 'FEMA_AE', 'FEMA_FSE', 'FEMA_FAE', 'FEMA_SAE', 'FEMA_FASE', 'fitlme_FE', 'fitlme_SE', 'fitlme_FSE', 'fitlme_par_FE', 'fitlme_par_SE', 'fitlme_par_FSE'};
results         = table('Size', [10 length(header)], 'VariableTypes', repmat({'double'}, length(header), 1), 'VariableNames', header);

% Loop over all number of observations and summarize
for vals = 1:length(settings.nObservations)
    toLoad = fullfile(rootDir, ['nObs-', num2str(settings.nObservations(vals), '%02d')]);
    
    % Load FEMA variables
    variables_FEMA_FE   = load(fullfile(toLoad, 'Results_FEMA_FE.mat'),     'FEMA_FE_elapsed');
    variables_FEMA_SE   = load(fullfile(toLoad, 'Results_FEMA_SE.mat'),     'FEMA_SE_elapsed');
    variables_FEMA_AE   = load(fullfile(toLoad, 'Results_FEMA_AE.mat'),     'FEMA_AE_elapsed');
    variables_FEMA_FSE  = load(fullfile(toLoad, 'Results_FEMA_FSE.mat'),    'FEMA_FSE_elapsed');
    variables_FEMA_FAE  = load(fullfile(toLoad, 'Results_FEMA_FAE.mat'),    'FEMA_FAE_elapsed');
    variables_FEMA_SAE  = load(fullfile(toLoad, 'Results_FEMA_SAE.mat'),    'FEMA_SAE_elapsed');
    variables_FEMA_FASE = load(fullfile(toLoad, 'Results_FEMA_FASE.mat'),   'FEMA_FASE_elapsed');
    
    % Load fitlme variables
    variables_lme_FE   = load(fullfile(toLoad, 'Results_fitlme_FE.mat'),  'lme_FE_Elapsed');
    variables_lme_SE   = load(fullfile(toLoad, 'Results_fitlme_SE.mat'),  'lme_SE_Elapsed');
    variables_lme_FSE  = load(fullfile(toLoad, 'Results_fitlme_FSE.mat'), 'lme_FSE_Elapsed');
    
    % Load parallel processing results
    variables_par_lme_FE = load(fullfile(toLoad, 'Results_fitlme_parTimes.mat'));
    
    % Save results
    results.NumObservations(vals)   = settings.nObservations(vals);
    results.FEMA_FE(vals)           = variables_FEMA_FE.FEMA_FE_elapsed;
    results.FEMA_SE(vals)           = variables_FEMA_SE.FEMA_SE_elapsed;
    results.FEMA_AE(vals)           = variables_FEMA_AE.FEMA_AE_elapsed;
    results.FEMA_FSE(vals)          = variables_FEMA_FSE.FEMA_FSE_elapsed;
    results.FEMA_FAE(vals)          = variables_FEMA_FAE.FEMA_FAE_elapsed;
    results.FEMA_SAE(vals)          = variables_FEMA_SAE.FEMA_SAE_elapsed;
    results.FEMA_FASE(vals)         = variables_FEMA_FASE.FEMA_FASE_elapsed;
    results.fitlme_FE(vals)         = variables_lme_FE.lme_FE_Elapsed;
    results.fitlme_SE(vals)         = variables_lme_SE.lme_SE_Elapsed;
    results.fitlme_FSE(vals)        = variables_lme_FSE.lme_FSE_Elapsed;
    results.fitlme_par_FE(vals)     = variables_par_lme_FE.lme_parFE_Elapsed;
    results.fitlme_par_SE(vals)     = variables_par_lme_FE.lme_parSE_Elapsed;
    results.fitlme_par_FSE(vals)    = variables_par_lme_FE.lme_parFSE_Elapsed;
end

% Save summary
save(fullfile(rootDir, 'Results_experiment03_nObs.mat'), 'results', 'settings');