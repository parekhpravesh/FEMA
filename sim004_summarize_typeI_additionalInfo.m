%% Summarize additional information from simulation 004: type I error rate
%% Initialize
rootDir              = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation004_Type1';
nyVars               = 50;
namesFEMA            = cellfun(@(x) strrep(x, ' ', ''), strcat({'FP_FEMA_y'}, num2str((1:nyVars)')), 'UniformOutput', false);
namesfitlmematrix    = cellfun(@(x) strrep(x, ' ', ''), strcat({'FP_fitlmematrix_y'}, num2str((1:nyVars)')), 'UniformOutput', false);
header               = [{'Repeat', 'FP_FEMA', 'FP_fitlmematrix'}, namesFEMA', namesfitlmematrix'];
locNamesFEMA         = not(cellfun(@isempty, regexpi(header, '^FP_FEMA_y')));
locNamesfitlmematrix = not(cellfun(@isempty, regexpi(header, '^FP_fitlmematrix_y')));

% Alpha values
alphaValues          = [0.05, 0.01, 1e-3, 1e-4];

% Initialize tables
results_005          = table('Size', [100 length(header)], 'VariableTypes', repmat({'double'}, length(header), 1), 'VariableNames', header);
results_001          = table('Size', [100 length(header)], 'VariableTypes', repmat({'double'}, length(header), 1), 'VariableNames', header);
results_0001         = table('Size', [100 length(header)], 'VariableTypes', repmat({'double'}, length(header), 1), 'VariableNames', header);
results_00001        = table('Size', [100 length(header)], 'VariableTypes', repmat({'double'}, length(header), 1), 'VariableNames', header);

% Some settings
nObs                 = 10000;
nRepeats             = 1:100;
nXvars               = 100;

%% For every repeat
for vals = 1:length(nRepeats)
    toLoad = fullfile(rootDir, ['RepNum-', num2str(nRepeats(vals),  '%05d'), '-nObs-', num2str(nObs, '%02d')]);
    
    % Load variables
    variables_FEMA_FSE  = load(fullfile(toLoad, 'Results_FEMA_FSE.mat'));
    variables_lme_FSE   = load(fullfile(toLoad, 'Results_fitlme_FSE.mat'));
    
    % Calculate p values and false positives from FEMA - get rid of first
    % column which is the contrast term
    pValues_FEMA     = (10.^(-sign(variables_FEMA_FSE.logpmat_FSE(2:end,:)) .* variables_FEMA_FSE.logpmat_FSE(2:end,:)))';
    altpValues_FEMA  = (normcdf(-abs(variables_FEMA_FSE.zmat_FSE(2:end,:)))*2)';
    alttpValues_FEMA = (2 * tcdf(-abs(variables_FEMA_FSE.zmat_FSE(2:end,:)), nObs-nXvars))';
    
    % p values from LME
    locPValues           = not(cellfun(@isempty, regexpi(variables_lme_FSE.results_FSE.Properties.VariableNames, '^p_')));
    pValues_fitlmematrix = variables_lme_FSE.results_FSE{:,locPValues};
    
    %% Alpha = 0.05
    alphaValue = alphaValues(1);
    
    % Number of false positives per y variable
    numFP_FEMA          = sum(altpValues_FEMA       < alphaValue, 2);
    numFP_fitlmematrix  = sum(pValues_fitlmematrix  < alphaValue, 2);
    
    % Overall number of false positives
    numFalsePositives_FEMA          = sum(altpValues_FEMA      < alphaValue, 'all');
    numFalsePositives_fitlmematrix  = sum(pValues_fitlmematrix < alphaValue, 'all');
    
    % Save results
    results_005.Repeat(vals)                = nRepeats(vals);
    results_005.FP_FEMA(vals)               = numFalsePositives_FEMA;
    results_005.FP_fitlmematrix(vals)       = numFalsePositives_fitlmematrix;
    results_005{vals, locNamesFEMA}         = numFP_FEMA';
    results_005{vals, locNamesfitlmematrix} = numFP_fitlmematrix';
    
    %% Alpha = 0.01
    alphaValue = alphaValues(2);
    
    % Number of false positives per y variable
    numFP_FEMA          = sum(altpValues_FEMA       < alphaValue, 2);
    numFP_fitlmematrix  = sum(pValues_fitlmematrix  < alphaValue, 2);
    
    % Overall number of false positives
    numFalsePositives_FEMA          = sum(altpValues_FEMA      < alphaValue, 'all');
    numFalsePositives_fitlmematrix  = sum(pValues_fitlmematrix < alphaValue, 'all');
    
    % Save results
    results_001.Repeat(vals)                = nRepeats(vals);
    results_001.FP_FEMA(vals)               = numFalsePositives_FEMA;
    results_001.FP_fitlmematrix(vals)       = numFalsePositives_fitlmematrix;
    results_001{vals, locNamesFEMA}         = numFP_FEMA';
    results_001{vals, locNamesfitlmematrix} = numFP_fitlmematrix';
    
    %% Alpha = 0.001
    alphaValue = alphaValues(3);
    
    % Number of false positives per y variable
    numFP_FEMA          = sum(altpValues_FEMA       < alphaValue, 2);
    numFP_fitlmematrix  = sum(pValues_fitlmematrix  < alphaValue, 2);
    
    % Overall number of false positives
    numFalsePositives_FEMA          = sum(altpValues_FEMA      < alphaValue, 'all');
    numFalsePositives_fitlmematrix  = sum(pValues_fitlmematrix < alphaValue, 'all');
    
    % Save results
    results_0001.Repeat(vals)                = nRepeats(vals);
    results_0001.FP_FEMA(vals)               = numFalsePositives_FEMA;
    results_0001.FP_fitlmematrix(vals)       = numFalsePositives_fitlmematrix;
    results_0001{vals, locNamesFEMA}         = numFP_FEMA';
    results_0001{vals, locNamesfitlmematrix} = numFP_fitlmematrix';
    
    %% Alpha = 0.0001
    alphaValue = alphaValues(4);
    
    % Number of false positives per y variable
    numFP_FEMA          = sum(altpValues_FEMA       < alphaValue, 2);
    numFP_fitlmematrix  = sum(pValues_fitlmematrix  < alphaValue, 2);
    
    % Overall number of false positives
    numFalsePositives_FEMA          = sum(altpValues_FEMA      < alphaValue, 'all');
    numFalsePositives_fitlmematrix  = sum(pValues_fitlmematrix < alphaValue, 'all');
    
    % Save results
    results_00001.Repeat(vals)                = nRepeats(vals);
    results_00001.FP_FEMA(vals)               = numFalsePositives_FEMA;
    results_00001.FP_fitlmematrix(vals)       = numFalsePositives_fitlmematrix;
    results_00001{vals, locNamesFEMA}         = numFP_FEMA';
    results_00001{vals, locNamesfitlmematrix} = numFP_fitlmematrix';
end

%% Average across repeats (of total number of false positives across y variables) rounded to positive inifnity
mean_total_FP_FEMA_005   = ceil(mean(results_005.FP_FEMA));
mean_total_FP_FEMA_001   = ceil(mean(results_001.FP_FEMA));
mean_total_FP_FEMA_0001  = ceil(mean(results_0001.FP_FEMA));
mean_total_FP_FEMA_00001 = ceil(mean(results_00001.FP_FEMA));

mean_total_FP_fitlmematrix_005   = ceil(mean(results_005.FP_fitlmematrix));
mean_total_FP_fitlmematrix_001   = ceil(mean(results_001.FP_fitlmematrix));
mean_total_FP_fitlmematrix_0001  = ceil(mean(results_0001.FP_fitlmematrix));
mean_total_FP_fitlmematrix_00001 = ceil(mean(results_00001.FP_fitlmematrix));

%% Save results
save(fullfile(rootDir, 'Results_experiment04_additionalInfo.mat'), 'results_*', 'mean_*');