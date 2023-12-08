%% Summarize results of simulation 001
% Set paths and get settings
rootDir  = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation001_ParamRecovery_sig2tvec1';
settings = load(fullfile(rootDir, 'Settings.mat'), 'settings');
settings = settings.settings;

% Initialize
colNames = {'Repeat', 'Bin', 'Voxel', 'TimeElapsed', ...
            'GTruth_X1', 'GTruth_X2', 'GTruth_X3', 'GTruth_X4', 'GTruth_X5', 'GTruth_F', 'GTruth_S', 'GTruth_E', 'GTruth_TotVar', ...
            'Est_X1',    'Est_X2',    'Est_X3',    'Est_X4',    'Est_X5',    'Est_F',    'Est_S',    'Est_E',    'Est_TotVar'};
count    = 1;

%% Loop over repeats and compile results
results = cell(settings.nRepeats, length(settings.allBins));
for repeats = 1:settings.nRepeats
    for bins = 1:length(settings.allBins)
        toLoad    = fullfile(rootDir, ['Repeat-', num2str(repeats, '%04d'), '-bin-', num2str(bins, '%04d')]);
        variables = load(fullfile(toLoad, 'Results.mat'), 'beta_hat', 'elapsed', 'sig2mat', 'sig2tvec');
        results{repeats, bins} = variables;
    end
end

%% Save results
save(fullfile(rootDir, 'Results_experiment01.mat'), 'results', 'settings');