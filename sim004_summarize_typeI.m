%% Summarize results of simulation 004: type I error rate

% Set paths and initialize
rootDir              = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/Simulation005_Type1';

% Various column names in summary
namesFEMA            = cellfun(@(x) strrep(x, ' ', ''), strcat({'FP_FEMA_y'}, num2str((1:20)')), 'UniformOutput', false);
namesfitlmematrix    = cellfun(@(x) strrep(x, ' ', ''), strcat({'FP_fitlmematrix_y'}, num2str((1:20)')), 'UniformOutput', false);
header               = [{'Repeat', 'FP_FEMA', 'FP_fitlmematrix'}, namesFEMA', namesfitlmematrix'];
locNamesFEMA         = not(cellfun(@isempty, regexpi(header, '^FP_FEMA_y')));
locNamesfitlmematrix = not(cellfun(@isempty, regexpi(header, '^FP_fitlmematrix_y')));

% Initialize
results              = table('Size', [100 length(header)], 'VariableTypes', ...
                             repmat({'double'}, length(header), 1), 'VariableNames', header);
nObs                 = 10000;
nRepeats             = 1:100;
alphaValue           = 0.05;
nXvars               = 100;

%% Loop over repeats and summarize
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
                        
    % Number of false positives per y variable
    numFP_FEMA          = sum(altpValues_FEMA       < alphaValue, 2);
    numFP_fitlmematrix  = sum(pValues_fitlmematrix  < alphaValue, 2);
    
    % Overall number of false positives
    numFalsePositives_FEMA          = sum(altpValues_FEMA      < alphaValue, 'all');
    numFalsePositives_fitlmematrix  = sum(pValues_fitlmematrix < alphaValue, 'all');
    
    % Save results
    results.Repeat(vals)                = nRepeats(vals);
    results.FP_FEMA(vals)               = numFalsePositives_FEMA;
    results.FP_fitlmematrix(vals)       = numFalsePositives_fitlmematrix;
    results{vals, locNamesFEMA}         = numFP_FEMA';
    results{vals, locNamesfitlmematrix} = numFP_fitlmematrix';
end

%% Save results
save('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/Simulation005_Type1/Results_experiment05.mat', 'results');