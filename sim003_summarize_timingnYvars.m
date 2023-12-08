%% Summarize simulation 003: timming as a function of number of y variables
% Set paths and initialize
rootDir         = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation003_CompareTimes_nYvars';
header          = {'NumYVars', 'FEMA_FSE', 'fitlme_FSE', 'fitlme_parFSE'};
results         = table('Size', [50 4], 'VariableTypes', repmat({'double'}, 4, 1), 'VariableNames', header);
nObs            = 10000;
nYVars          = 100:100:5000;

% Loop over all number of y variables and summarize
for vals = 1:length(nYVars)
    toLoad = fullfile(rootDir, ['nyVars-', num2str(nYVars(vals),  '%05d'), '-nObs-', num2str(nObs, '%02d')]);
    
    % Load variables
    variables_FEMA_FSE  = load(fullfile(toLoad, 'Results_FEMA_FSE.mat'),      'FEMA_FSE_elapsed');
    variables_lme_FSE   = load(fullfile(toLoad, 'Results_fitlme_FSE.mat'),    'lme_FSE_Elapsed');
    variables_lme_par   = load(fullfile(toLoad, 'Results_fitlme_parFSE.mat'), 'lme_parFSE_Elapsed');
    
    % Save results
    results.NumYVars(vals)      = nYVars(vals);
    results.FEMA_FSE(vals)      = variables_FEMA_FSE.FEMA_FSE_elapsed;
    results.fitlme_FSE(vals)    = variables_lme_FSE.lme_FSE_Elapsed;
    results.fitlme_parFSE(vals) = variables_lme_par.lme_parFSE_Elapsed;
end

% Save summary
save(fullfile(rootDir, 'Results_experiment03_nYvars.mat'), 'results');