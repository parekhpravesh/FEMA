%% Simulation 001: effect of binning and parameter recovery
% In this experiment, we will generate 50 repeats of 2000 y variables and
% 10,000 observations. There will be 5 X variables and the ground truth 
% effect sizes will be the same across all repeats. 
% For each repeat, we will generate new X variables and new y variables. 
% Parameter recovery will be using method of moments (MoM) estimator from
% FEMA and will be for 20 different values of bins

%% Include relevant paths
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/codebase/cmig_tools-2.3.0/cmig_tools_utils/matlab');
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/codebase/cmig_tools-2.3.0/FEMA');

%% Prepare output directory
outDir = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation001_ParamRecovery_sig2tvec1';
if not(exist(outDir, 'dir'))
    mkdir(outDir);
end

%% Setup parameters 
nObservations    = 10000;
nXvars           = 5;
nyVars           = 2000;
nFamMembers      = 5;
nRepObservations = 5;
RandomEffects    = {'F', 'S', 'E'}; 
eid              = ones(nObservations, 1);
agevec           = zeros(nObservations, 1);
y_RFX            = nan(nObservations, nyVars);
epsmin           = 0.2; 
sigmaLow         = 0.2;
nRepeats         = 50;
allBins          = [1:5, 10:10:50, unique(round(logspace(2,log10(nyVars),10),0))];
niter            = 1;

%% Control random number here and generate ground truth effect sizes
rng(20230202, 'twister');

% Beta is uniformally distributed in the range -0.02 and 0.02 => rand
betaLow     = -0.2;
betaHigh    = 0.2;
beta        = betaLow + (betaHigh-betaLow).*rand(nXvars,nyVars);

% Random effects are uniformally distributed in the range 0.2 and 0.8
% Note that the range is the sum of F and S and does not imply that each of
% the F and S values will be at least equal to 0.2
sig2mat_true = NaN(length(RandomEffects)-1, size(y_RFX,2)); 
ivec         = 1:size(y_RFX,2);
while ~isempty(ivec)
    sig2mat_true(:,ivec) = rand([size(sig2mat_true,1) length(ivec)]);
    cond                 = (sum(sig2mat_true,1) >= sigmaLow) & (sum(sig2mat_true,1) < 1-epsmin);
    ivec                 = find(not(cond));
    logging('length(ivec)=%d',length(ivec));
    if isempty(ivec)
        break;
    end
end
sig2mat_true = cat(1, sig2mat_true, 1-sum(sig2mat_true,1));

% Total residual variance is 1
sig2tvec_true = 1;

%% Generate seeds
allSeeds = randi(999999, nRepeats, 1);

%% Summary of all ground truth and settings
settings.nObservations          = nObservations;
settings.nXvars                 = nXvars;
settings.nyVars                 = nyVars;
settings.nFamMembers            = nFamMembers;
settings.nRepObservations       = nRepObservations;
settings.RandomEffects          = RandomEffects; 
settings.epsmin                 = epsmin; 
settings.sigmaLow               = sigmaLow;
settings.nRepeats               = nRepeats;
settings.allBins                = allBins;
settings.niter                  = niter;
settings.allSeeds               = allSeeds;
settings.GTruth.betaLow         = betaLow;
settings.GTruth.betaHigh        = betaHigh;
settings.GTruth.beta            = beta;
settings.GTruth.sig2mat_true    = sig2mat_true;
settings.GTruth.sig2tvec_true   = sig2tvec_true;

save(fullfile(outDir, 'Settings.mat'), 'settings');

%% Start main repeat loop
for repeats = 1:nRepeats
    
    % Set seed for this repeat
    rng(allSeeds(repeats), 'twister');

    % Generate X variables: X is normally distributed => randn
    X     = randn(nObservations, nXvars);
    y_FFX = X * beta;

    % Genereate IDs - control number of members in family and number of 
    % % up to 5 observations per individual
    iid_int = sort(randsample(repmat(1:nObservations, 1, nRepObservations), nObservations));

    % up to 5 individuals per family
    fid_int = ceil(iid_int / nFamMembers); 

    % Initialize
    iid = cell(size(iid_int)); 
    fid = cell(size(fid_int)); 

    % Generate individual and family IDs
    for i=1:length(iid_int)
        iid{i}=sprintf('I%i', iid_int(i)); 
    end
    for i=1:length(fid_int)
        fid{i}=sprintf('F%i', fid_int(i)); 
    end

    % Using FSE model, generate y_RFX
    % Parse family structure
    clusterinfo = FEMA_parse_family(iid, eid, fid, agevec, [], ...
                                   'RandomEffects', RandomEffects);

    % Generate y_RFX
    nfam = length(clusterinfo); 
    for fi = 1:nfam
        if mod(fi,100)==0
            logging('fi=%d/%d',fi,nfam);
        end
        tmp = 0;
        for ri = 1:length(RandomEffects)
            tmp = tmp + sqrt(sig2mat_true(ri,:)) .*                        ...
                        mvnrnd(zeros(length(clusterinfo{fi}.jvec_fam), 1), ...
                        double(getfield(clusterinfo{fi},                   ...
                        sprintf('V_%s',RandomEffects{ri}))),nyVars)'; %#ok<GFLD>
        end
        y_RFX(clusterinfo{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
    end

    % Put together as a single y variable
    ymat = y_FFX + y_RFX;
    
    % Start loop for every bin size
    for bins = 1:length(allBins)
        
        % Skip if results exist
        tmpOutDir = fullfile(outDir, ['Repeat-', num2str(repeats, '%04d'), '-bin-', num2str(bins, '%04d')]);
        if not(exist(tmpOutDir, 'dir'))
            mkdir(tmpOutDir);
        end
        if exist(fullfile(tmpOutDir, 'Results.mat'), 'file')
            warning(['Skipping: Repeat ', num2str(repeats, '%04d'), ', bin ', num2str(bins, '%04d')]);
            continue;
        end
        
        % Initialize timer
        timeInit = tic;

        % Estimate model parameters using FEMA
        [beta_hat,      beta_se,        zmat,        logpmat,              ...
         sig2tvec,      sig2mat,        ~,           ~,                    ...
         ~,             ~,              ~,           ~,                    ...
         ~,             ~,              binvec_save, ~,                    ...
         ~,             ~,              ~] =                               ...
        FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),    ...
                allBins(bins), [], 'RandomEffects', RandomEffects, 'returnReusable', true);
            
        % End timer
        elapsed = toc(timeInit);
        
        % Save everything
        save(fullfile(tmpOutDir, 'Results.mat'), 'beta_hat', 'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat', ...
                                                 'elapsed', 'repeats', 'bins', 'ymat', 'X', 'iid', 'fid', 'binvec_save');
                            
        % Update user
        disp(['Finished: Repeat ', num2str(repeats, '%04d'), ', bin ', num2str(bins, '%04d'), ' in ', num2str(toc(timeInit), '%.2f') 's']);
    end
end