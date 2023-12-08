%% Summarize results of simulation 006: type I error rate as a function of bins
% Set paths and initialize
rootDir = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation006_Type1Bins';

% Some settings
nObs            = 10000;
nRepeats        = 1:100;
alphaValue      = 0.05;
nXvars          = 100;
allBins         = 1:30;

% Various column names in summary
nyVars          = 500;
namesFEMA       = cellfun(@(x) strrep(x, ' ', ''), strcat({'FP_FEMA_y'}, num2str((1:nyVars)')), 'UniformOutput', false);
namesRepeats    = cellfun(@(x) strrep(x, ' ', ''), strcat({'Repeat'},    num2str((nRepeats)')), 'UniformOutput', false);
header          = [{'BinSize', 'Repeat', 'FP_FEMA'}, namesFEMA'];
header_uqBins   = {'BinSize', 'Repeat', 'NumTotalBins'};
header_maxY     = {'BinSize', 'Repeat', 'MaxYVars'};
locNamesFEMA    = not(cellfun(@isempty, regexpi(header, '^FP_FEMA_y')));

% Initialize
results         = table('Size', [length(nRepeats)*length(allBins) length(header)], 'VariableTypes', ...
                        repmat({'double'}, length(header), 1), 'VariableNames', header);
uqBins          = table('Size', [length(nRepeats)*length(allBins) length(header_uqBins)], 'VariableTypes', ...
                        repmat({'double'}, length(header_uqBins), 1), 'VariableNames', header_uqBins);
maxYvars        = table('Size', [length(nRepeats)*length(allBins) length(header_maxY)], 'VariableTypes', ...
                        repmat({'double'}, length(header_maxY), 1), 'VariableNames', header_maxY);

%% Loop over bins and repeats and summarize
count = 1;
for bins = 1:length(allBins)
    currBin = allBins(bins);
    for repeats = 1:length(nRepeats)
        currDir = fullfile(rootDir, ['RepNum-', num2str(repeats, '%05d'), '-nObs-', num2str(nObs, '%04d')]);
        toLoad  = fullfile(currDir, ['Results_FEMA_FSE_bin-', num2str(currBin, '%02d'), '.mat']);
        
        % Load variables
        variables_FEMA_FSE = load(toLoad, 'zmat_FSE', 'logpmat_FSE', 'binvec_save');
        
        % Calculate p values and false positives
        pValues_FEMA = (10.^(-sign(variables_FEMA_FSE.zmat_FSE) .* variables_FEMA_FSE.logpmat_FSE))';
        % pValues_FEMA     = (10.^(-sign(variables_FEMA_FSE.logpmat_FSE) .* variables_FEMA_FSE.logpmat_FSE))';
        % altpValues_FEMA  = (normcdf(-abs(variables_FEMA_FSE.zmat_FSE))*2)';
        % alttpValues_FEMA = (2 * tcdf(-abs(variables_FEMA_FSE.zmat_FSE), nObs-nXvars))';
                                    
        % Number of false positives per y variable
        numFP_FEMA = sum(pValues_FEMA < alphaValue, 2);
        
        % Overall number of false positives
        numFalsePositives_FEMA = sum(pValues_FEMA < alphaValue, 'all');
        
        % Save results
        results.BinSize(count)       = currBin;
        results.Repeat(count)        = nRepeats(repeats);
        results.FP_FEMA(count)       = numFalsePositives_FEMA;
        results{count, locNamesFEMA} = numFP_FEMA';

        % Count the number of unique bins
        uqBins.BinSize(count)        = currBin;
        uqBins.Repeat(count)         = nRepeats(repeats);
        uqBins.NumTotalBins(count)   = length(unique(variables_FEMA_FSE.binvec_save));

        % Find the max number of y variables assigned to any bin
        [a, b]                       = histcounts(categorical(variables_FEMA_FSE.binvec_save));
        maxYvars.BinSize(count)      = currBin;
        maxYvars.Repeat(count)       = nRepeats(repeats);
        maxYvars.MaxYVars(count)     = max(a);

        % Update counter
        count = count + 1;

        % Clear up
        clear toLoad variables_FEMA_FSE pValues_FEMA
    end

    % Update
    disp(['Finished bin ', num2str(currBin)]);
end

%% Save results
save(fullfile(rootDir, 'Results_experiment06.mat'), 'results', 'maxYvars', 'uqBins');