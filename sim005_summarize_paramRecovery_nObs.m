%% Summarize results of simulation 005: parameter recovery as a function of number of observations
% Set paths and initialize
rootDir         = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation005_Compare_fitlme_perms_nn_nObs';
settings        = load(fullfile(rootDir, 'Settings.mat'), 'settings');
settings        = settings.settings;
numRepeats      = 5;
numRFX          = length(settings.RandomEffects);
est_FFX_FEMA    = zeros(settings.nXvars, settings.nyVars, numRepeats, length(settings.nObservations));
est_FFX_LME     = zeros(settings.nXvars, settings.nyVars, numRepeats, length(settings.nObservations));
est_RFX_FEMA    = zeros(numRFX,          settings.nyVars, numRepeats, length(settings.nObservations));
est_RFX_LME     = zeros(numRFX,          settings.nyVars, numRepeats, length(settings.nObservations));

%% Loop over every number of observation and summarize
for vals = 1:length(settings.nObservations)
    for repeats     = 1:numRepeats
        toLoad      = fullfile(rootDir, [num2str(vals, '%04d'), '-nObs-', num2str(settings.nObservations(vals), '%04d')], [num2str(vals, '%04d'), '-nObs-', num2str(settings.nObservations(vals), '%04d'), '-Repeat-', num2str(repeats, '%02d')]);
        vars_FEMA   = load(fullfile(toLoad, 'Results_FEMA.mat'), 'beta_*', 'FEMAelapsed', 'sig2*');
        vars_fitlme = load(fullfile(toLoad, 'Results_fitlmematrix_summary.mat'), 'variables_lme');
        est_FFX_FEMA(1:settings.nXvars, 1:settings.nyVars, repeats, vals) = vars_FEMA.beta_hat(2:end,:);
        est_FFX_LME(1:settings.nXvars,  1:settings.nyVars, repeats, vals) = vars_fitlme.variables_lme.beta_hat_fitlme;
        est_RFX_FEMA(1:numRFX,          1:settings.nyVars, repeats, vals) = vars_FEMA.sig2mat;
        est_RFX_LME(1:numRFX,           1:settings.nyVars, repeats, vals) = vars_fitlme.variables_lme.sig2mat_fitlme.^2;
    end
end

%% Sum(Sum(Average(Squared difference between ground truth and estimated value, over repeats), over y variables), over X variables)
TotalmeanSqDiff_GTruth_FEMA_FFX = sum(squeeze(sum(squeeze(mean((settings.GTruth.beta - est_FFX_FEMA).^2, 3)),2)), 1);
TotalmeanSqDiff_GTruth_LME_FFX  = sum(squeeze(sum(squeeze(mean((settings.GTruth.beta - est_FFX_LME).^2,  3)),2)), 1);

TotalmeanSqDiff_GTruth_FEMA_RFX = sum(squeeze(sum(squeeze(mean((settings.GTruth.sig2mat_true - est_RFX_FEMA).^2, 3)),2)), 1);
TotalmeanSqDiff_GTruth_LME_RFX  = sum(squeeze(sum(squeeze(mean((settings.GTruth.sig2mat_true - est_RFX_LME).^2,  3)),2)), 1);

%% Sum(Sum(Squared difference between ground truth and estimated value, over y variables), over X variables)
TotalSqDiff_GTruth_FEMA_FFX = squeeze(sum(sum((settings.GTruth.beta - est_FFX_FEMA).^2, 2), 1));
TotalSqDiff_GTruth_LME_FFX  = squeeze(sum(sum((settings.GTruth.beta - est_FFX_LME).^2,  2), 1));

TotalSqDiff_GTruth_FEMA_RFX = squeeze(sum(sum((settings.GTruth.sig2mat_true - est_RFX_FEMA).^2, 2), 1));
TotalSqDiff_GTruth_LME_RFX  = squeeze(sum(sum((settings.GTruth.sig2mat_true - est_RFX_LME).^2,  2), 1));

%% sum(Squared difference between ground truth and estimated value, over y variables)
TotalySqDiff_GTruth_FEMA_FFX = squeeze(sum((settings.GTruth.beta - est_FFX_FEMA).^2, 2));
TotalySqDiff_GTruth_LME_FFX  = squeeze(sum((settings.GTruth.beta - est_FFX_LME).^2, 2));

TotalySqDiff_GTruth_FEMA_RFX = squeeze(sum((settings.GTruth.sig2mat_true - est_RFX_FEMA).^2, 2));
TotalySqDiff_GTruth_LME_RFX  = squeeze(sum((settings.GTruth.sig2mat_true - est_RFX_LME).^2, 2));

%% Save results
save(fullfile(rootDir, 'Results_experiment05.mat'));