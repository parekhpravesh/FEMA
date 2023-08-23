%% Simulation 003: comparison with fitlme - time focused (function of number of y variables)
%% Include relevant paths
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/cmig_tools_internal-beta/cmig_tools_utils/matlab');
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/cmig_tools_internal-beta/FEMA');

%% Prepare output directory
outDir = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/Simulation004_Compare_times_yVars_fitlmematrix';
if not(exist(outDir, 'dir'))
    mkdir(outDir);
end

%% Setup parameters 
nObservations    = 10000;
nXvars           = 5;
nyVars           = 100:100:5000;
nFamMembers      = 5;
nRepObservations = 5;
RandomEffects    = {'F', 'S', 'E'}; 
epsmin           = 0.2; 
sigmaLow         = 0.2;
allBins          = 20;
niter            = 1;

%% Generate seeds
rng(20230216, 'twister');
allSeeds = randi(99999999, length(nyVars), 1);

%% Header
header = {'YvarNum', 'Beta_X1', 'Beta_X2', 'Beta_X3', 'Beta_X4', 'Beta_X5', ...
                     'SE_X1',   'SE_X2',   'SE_X3',   'SE_X4',   'SE_X5',   ...
                     'V_F',     'V_S',     'V_E',                           ...
                     'Low_VF',  'Low_VS',  'Low_VE',                        ...
                     'Upp_VF',  'Upp_VS',  'Upp_VE'};

%% Initialize parallel computing
% local            = parcluster('local');
% local.NumThreads = 2;
% pool             = local.parpool(20, 'IdleTimeout', 240);

%% Loop over number of y variables
for loopVar = 1:length(nyVars)

    % Control random number here and generate ground truth effect sizes
    rng(allSeeds(loopVar), 'twister');
    
    % Current number of y variables
    currNyVars = nyVars(loopVar);

    % Beta is uniformally distributed in the range -0.02 and 0.02 => rand
    betaLow     = -0.2;
    betaHigh    = 0.2;
    beta        = betaLow + (betaHigh-betaLow).*rand(nXvars,currNyVars);

    % Random effects are uniformally distributed in the range 0.2 and 0.8
    % Note that the range is the sum of F and S and does not imply that each of
    % the F and S values will be at least equal to 0.2
    sig2mat_true = NaN(length(RandomEffects)-1, currNyVars); 
    ivec         = 1:currNyVars;
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

    % Summary of all ground truth and settings
    settings.nObservations          = nObservations;
    settings.nXvars                 = nXvars;
    settings.nyVars                 = currNyVars;
    settings.nFamMembers            = nFamMembers;
    settings.nRepObservations       = nRepObservations;
    settings.RandomEffects          = RandomEffects; 
    settings.epsmin                 = epsmin; 
    settings.sigmaLow               = sigmaLow;
    settings.allBins                = allBins;
    settings.niter                  = niter;
    settings.allSeeds               = allSeeds;
    settings.GTruth.betaLow         = betaLow;
    settings.GTruth.betaHigh        = betaHigh;
    settings.GTruth.beta            = beta;
    settings.GTruth.sig2mat_true    = sig2mat_true;
    settings.GTruth.sig2tvec_true   = sig2tvec_true;

    save(fullfile(outDir, ['Settings_numYvars_', num2str(currNyVars, '%05d') '.mat']), 'settings');
    
    % Current number of observations
    currObs = nObservations(1);
    
    % Initialize a few variables
    eid              = ones(currObs, 1);
    agevec           = zeros(currObs, 1);
    y_RFX            = nan(currObs, currNyVars);
    
    % Generate X variables: X is normally distributed => randn
    % Make the first X variable as the intercept
    X     = [ones(currObs,1), randn(currObs, nXvars-1)];
    y_FFX = X * beta;

    % Genereate IDs - control number of members in family and number of
    % % up to 5 observations per individual
    iid_int = sort(randsample(repmat(1:currObs, 1, nRepObservations), currObs));

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
                        sprintf('V_%s',RandomEffects{ri}))),currNyVars)'; %#ok<GFLD>
        end
        y_RFX(clusterinfo{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
    end

    % Put together as a single y variable
    ymat = y_FFX + y_RFX;
    
    % Output directory
    tmpOutDir = fullfile(outDir, ['nyVars-', num2str(currNyVars, '%05d'), '-nObs-', num2str(currObs, '%04d')]);
    if not(exist(tmpOutDir, 'dir'))
        mkdir(tmpOutDir);
    end
    for i=1:length(fid_int)
        fid{i}=sprintf('F%i', fid_int(i));
    end
    
    % Fit FSE models
    % Initialize timer for FEMA
    FEMA_FSE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_FSE, beta_se_FSE, zmat_FSE, logpmat_FSE, sig2tvec_FSE, sig2mat_FSE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),                ...
              allBins, [], 'RandomEffects', {'F', 'S', 'E'});
            
    % End timer for FEMA
    FEMA_FSE_elapsed = toc(FEMA_FSE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_FSE.mat'), 'beta_hat_FSE',    'beta_se_FSE', 'zmat_FSE', 'logpmat_FSE', 'sig2tvec_FSE', 'sig2mat_FSE', ...
                                                      'FEMA_FSE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                         
    % Now loop over all phenotypes and estimate model using fitlmematrix
    % Ignore the first X column as intercept
    mdls_FSE      = cell(currNyVars, 1);
    lme_FSE_Init  = tic;
    for yVar = 1:currNyVars
        mdls_FSE{yVar} = fitlmematrix(X, ymat(:,yVar), {ones(length(fid),1), ones(length(iid),1)}, {fid', iid'});
    end
    lme_FSE_Elapsed = toc(lme_FSE_Init);
    
    % Extract out everything relevant
    results_FSE = table('Size', [currNyVars, 20], 'VariableTypes', repmat({'double'}, 20, 1), 'VariableNames', header);
    for yVar = 1:currNyVars
        [a, b, c]                   = covarianceParameters(mdls_FSE{yVar});
        results_FSE.YvarNum(yVar)   = yVar;
        results_FSE.Beta_X1(yVar)   = mdls_FSE{yVar}.Coefficients.Estimate(1);
        results_FSE.Beta_X2(yVar)   = mdls_FSE{yVar}.Coefficients.Estimate(2);
        results_FSE.Beta_X3(yVar)   = mdls_FSE{yVar}.Coefficients.Estimate(3);
        results_FSE.Beta_X4(yVar)   = mdls_FSE{yVar}.Coefficients.Estimate(4);
        results_FSE.Beta_X5(yVar)   = mdls_FSE{yVar}.Coefficients.Estimate(5);
        results_FSE.SE_X1(yVar)     = mdls_FSE{yVar}.Coefficients.SE(1);
        results_FSE.SE_X2(yVar)     = mdls_FSE{yVar}.Coefficients.SE(2);
        results_FSE.SE_X3(yVar)     = mdls_FSE{yVar}.Coefficients.SE(3);
        results_FSE.SE_X4(yVar)     = mdls_FSE{yVar}.Coefficients.SE(4);
        results_FSE.SE_X5(yVar)     = mdls_FSE{yVar}.Coefficients.SE(5);
        results_FSE.V_F(yVar)       = c{1}.Estimate;
        results_FSE.V_S(yVar)       = c{2}.Estimate;
        results_FSE.V_E(yVar)       = c{3}.Estimate;
        results_FSE.Low_VF(yVar)    = c{1}.Lower;
        results_FSE.Low_VS(yVar)    = c{2}.Lower;
        results_FSE.Low_VE(yVar)    = c{3}.Lower;
        results_FSE.Upp_VF(yVar)    = c{1}.Upper;
        results_FSE.Upp_VS(yVar)    = c{2}.Upper;
        results_FSE.Upp_VE(yVar)    = c{3}.Upper;
    end
    
    % Save everything
    save(fullfile(tmpOutDir, 'Results_fitlme_FSE.mat'), 'results_FSE', 'lme_FSE_Elapsed');
    
    % Start a pool
    local            = parcluster('local');
    local.NumThreads = 2;
    pool             = local.parpool(20, 'IdleTimeout', 240);
    
    % Redo processing with parallel computing
    mdls_parFSE      = cell(currNyVars, 1);
    lme_parFSE_Init  = tic;
    Zvariable        = {ones(length(fid),1), ones(length(iid),1)};
    Gvariable        = {fid', iid'};
    parfor yVar = 1:currNyVars
        mdls_parFSE{yVar} = fitlmematrix(X, ymat(:,yVar), Zvariable, Gvariable);
    end
    lme_parFSE_Elapsed = toc(lme_parFSE_Init);
    
    % Delete the pool
    delete(pool);
    
    % Save timing
    save(fullfile(tmpOutDir, 'Results_fitlme_parFSE.mat'), 'lme_parFSE_Elapsed');
    
%     % Extract out everything relevant
%     results_parFSE = table('Size', [currNyVars, 20], 'VariableTypes', repmat({'double'}, 20, 1), 'VariableNames', header);
%     for yVar = 1:currNyVars
%         [a, b, c]                      = covarianceParameters(mdls_parFSE{yVar});
%         results_parFSE.YvarNum(yVar)   = yVar;
%         results_parFSE.Beta_X1(yVar)   = mdls_parFSE{yVar}.Coefficients.Estimate(1);
%         results_parFSE.Beta_X2(yVar)   = mdls_parFSE{yVar}.Coefficients.Estimate(2);
%         results_parFSE.Beta_X3(yVar)   = mdls_parFSE{yVar}.Coefficients.Estimate(3);
%         results_parFSE.Beta_X4(yVar)   = mdls_parFSE{yVar}.Coefficients.Estimate(4);
%         results_parFSE.Beta_X5(yVar)   = mdls_parFSE{yVar}.Coefficients.Estimate(5);
%         results_parFSE.SE_X1(yVar)     = mdls_parFSE{yVar}.Coefficients.SE(1);
%         results_parFSE.SE_X2(yVar)     = mdls_parFSE{yVar}.Coefficients.SE(2);
%         results_parFSE.SE_X3(yVar)     = mdls_parFSE{yVar}.Coefficients.SE(3);
%         results_parFSE.SE_X4(yVar)     = mdls_parFSE{yVar}.Coefficients.SE(4);
%         results_parFSE.SE_X5(yVar)     = mdls_parFSE{yVar}.Coefficients.SE(5);
%         results_parFSE.V_F(yVar)       = c{1}.Estimate;
%         results_parFSE.V_S(yVar)       = c{2}.Estimate;
%         results_parFSE.V_E(yVar)       = c{3}.Estimate;
%         results_parFSE.Low_VF(yVar)    = c{1}.Lower;
%         results_parFSE.Low_VS(yVar)    = c{2}.Lower;
%         results_parFSE.Low_VE(yVar)    = c{3}.Lower;
%         results_parFSE.Upp_VF(yVar)    = c{1}.Upper;
%         results_parFSE.Upp_VS(yVar)    = c{2}.Upper;
%         results_parFSE.Upp_VE(yVar)    = c{3}.Upper;
%     end
%     
    % Update user
    disp(['Finished: Parallel FSE: Number of y variables ', num2str(currNyVars, '%04d')]);
end