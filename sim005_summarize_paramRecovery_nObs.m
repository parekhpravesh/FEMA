%% Summarize results of simulation 005: parameter recovery as a function of number of observations

% Set paths and initialize
rootDir                 = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-08-02_Redone/Simulation004_Compare_fitlme_perms_nn_repeats';
settings                = load(fullfile(rootDir, 'Settings.mat'), 'settings');
settings                = settings.settings;
SqDiff_FEMA_truth       = zeros(5, 500, 16);
SqDiff_LME_truth        = zeros(5, 500, 16);
SqDiff_FEMA_LME         = zeros(5, 500, 16);
SqDiff_RFX_FEMA_truth   = zeros(3, 500, 16);
SqDiff_RFX_LME_truth    = zeros(3, 500, 16);
SqDiff_RFX_FEMA_LME     = zeros(3, 500, 16);

%% Loop over every number of observation and summarize
for vals = 1:16
    tmp_SqDiff_FEMA_truth = zeros(5, 500, 5);
    tmp_SqDiff_LME_truth  = zeros(5, 500, 5);
    tmp_SqDiff_FEMA_LME   = zeros(5, 500, 5);
    
    tmp_RFX_SqDiff_FEMA_truth = zeros(3, 500, 5);
    tmp_RFX_SqDiff_LME_truth  = zeros(3, 500, 5);
    tmp_RFX_SqDiff_FEMA_LME   = zeros(3, 500, 5);

    for repeats     = 1:5
        toLoad      = fullfile(rootDir, [num2str(vals, '%04d'), '-nObs-', num2str(settings.nObservations(vals), '%04d')], [num2str(vals, '%04d'), '-nObs-', num2str(settings.nObservations(vals), '%04d'), '-Repeat-', num2str(repeats, '%02d')]);
        vars_FEMA   = load(fullfile(toLoad, 'Results_FEMA.mat'), 'beta_*', 'FEMAelapsed', 'sig2*');
        vars_fitlme = load(fullfile(toLoad, 'Results_fitlmematrix_summary.mat'), 'variables_lme');
        
        % Squared difference
        tmp_SqDiff_FEMA_truth(1:5, :, repeats) = (settings.GTruth.beta        - vars_FEMA.beta_hat(2:end,:)).^2;
        tmp_SqDiff_LME_truth(1:5,  :, repeats) = (settings.GTruth.beta        - vars_fitlme.variables_lme.beta_hat_fitlme).^2;
        tmp_SqDiff_FEMA_LME(1:5,   :, repeats) = (vars_FEMA.beta_hat(2:end,:) - vars_fitlme.variables_lme.beta_hat_fitlme).^2;
        
        % Squared difference - random effects; remember to square
        % fitlmematrix estimates before comparing
        tmp_RFX_SqDiff_FEMA_truth(1:3, :, repeats) = (settings.GTruth.sig2mat_true - vars_FEMA.sig2mat).^2;
        tmp_RFX_SqDiff_LME_truth(1:3, :, repeats)  = (settings.GTruth.sig2mat_true - (vars_fitlme.variables_lme.sig2mat_fitlme).^2).^2;
        tmp_RFX_SqDiff_FEMA_LME(1:3, :, repeats)   = (vars_FEMA.sig2mat            - (vars_fitlme.variables_lme.sig2mat_fitlme).^2).^2;
    end
    
    % Average squared difference across repeats
    SqDiff_FEMA_truth(:, :, vals) = mean(tmp_SqDiff_FEMA_truth, 3);
    SqDiff_LME_truth(:, :, vals)  = mean(tmp_SqDiff_LME_truth,  3);
    SqDiff_FEMA_LME(:, :, vals)   = mean(tmp_SqDiff_FEMA_LME,   3);
    
    % Average squared difference across repeats - random effects
    SqDiff_RFX_FEMA_truth(:, :, vals) = mean(tmp_RFX_SqDiff_FEMA_truth, 3);
    SqDiff_RFX_LME_truth(:,  :, vals) = mean(tmp_RFX_SqDiff_LME_truth,  3);
    SqDiff_RFX_FEMA_LME(:,   :, vals) = mean(tmp_RFX_SqDiff_FEMA_LME,   3);
    
end

%% Total average squared differences
totalSqDiff_FEMA_truth = squeeze(sum(SqDiff_FEMA_truth,2));
totalSqDiff_LME_truth  = squeeze(sum(SqDiff_LME_truth, 2));
totalSqDiff_FEMA_LME   = squeeze(sum(SqDiff_FEMA_LME,  2));

% Total average squared difference - random effects
totalSqDiff_RFX_FEMA_truth  = squeeze(sum(SqDiff_RFX_FEMA_truth, 2));
totalSqDiff_RFX_LME_truth   = squeeze(sum(SqDiff_RFX_LME_truth,  2));
totalSqDiff_RFX_FEMA_LME    = squeeze(sum(SqDiff_RFX_FEMA_LME,   2));

%% Average of the total squared difference across x variables
AvgTotalSqDiff_FEMA_truth = mean(totalSqDiff_FEMA_truth, 1);
AvgTotalSqDiff_LME_truth  = mean(totalSqDiff_LME_truth,  1);
AvgTotalSqDiff_FEMA_LME   = mean(totalSqDiff_FEMA_LME,   1);

% Average of the total squared difference across x variables - random
% effects
AvgTotalSqDiff_RFX_FEMA_truth   = mean(totalSqDiff_RFX_FEMA_truth, 1);
AvgTotalSqDiff_RFX_LME_truth    = mean(totalSqDiff_RFX_LME_truth,  1);
AvgTotalSqDiff_RFX_FEMA_LME     = mean(totalSqDiff_RFX_FEMA_LME,   1);

%% Sum of the total squared difference across x variables
SumTotalSqDiff_FEMA_truth = sum(totalSqDiff_FEMA_truth, 1);
SumTotalSqDiff_LME_truth  = sum(totalSqDiff_LME_truth,  1);
SumTotalSqDiff_FEMA_LME   = sum(totalSqDiff_FEMA_LME,   1);

% Sum of the total squared difference across x variables - random effects
SumTotalSqDiff_RFX_FEMA_truth = sum(totalSqDiff_RFX_FEMA_truth, 1);
SumTotalSqDiff_RFX_LME_truth  = sum(totalSqDiff_RFX_LME_truth,  1);
SumTotalSqDiff_RFX_FEMA_LME   = sum(totalSqDiff_RFX_FEMA_LME,   1);

%% Save results
save('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-08-02_Redone/Simulation004_Compare_fitlme_perms_nn_repeats/Summary_Simulation004_Compare_fitlme_perms_nn_repeats.mat');