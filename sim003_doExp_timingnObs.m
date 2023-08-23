%% Simulation 003: comparison with fitlme - time focused (function of number of observations)
%% Include relevant paths
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/cmig_tools_internal-beta/cmig_tools_utils/matlab');
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/cmig_tools_internal-beta/FEMA');

%% Prepare output directory
outDir = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/Simulation003_Compare_times_again';
if not(exist(outDir, 'dir'))
    mkdir(outDir);
end

%% Setup parameters 
nObservations    = linspace(2000, 20000, 10);
nXvars           = 5;
nyVars           = 50;
nFamMembers      = 5;
nRepObservations = 5;
RandomEffects    = {'F', 'A', 'S', 'E'}; 
epsmin           = 0.2; 
sigmaLow         = 0.2;
allBins          = 20;
niter            = 1;

%% Header
header = {'YvarNum', 'Beta_X1', 'Beta_X2', 'Beta_X3', 'Beta_X4', 'Beta_X5', ...
                     'SE_X1',   'SE_X2',   'SE_X3',   'SE_X4',   'SE_X5',   ...
                     'V_F',     'V_S',     'V_E',                           ...
                     'Low_VF',  'Low_VS',  'Low_VE',                        ...
                     'Upp_VF',  'Upp_VS',  'Upp_VE'};
                 
%% Control random number here and generate ground truth effect sizes
rng(20230207, 'twister');

% Beta is uniformally distributed in the range -0.02 and 0.02 => rand
betaLow     = -0.2;
betaHigh    = 0.2;
beta        = betaLow + (betaHigh-betaLow).*rand(nXvars,nyVars);

% Random effects are uniformally distributed in the range 0.2 and 0.8
% Note that the range is the sum of F and S and does not imply that each of
% the F and S values will be at least equal to 0.2
sig2mat_true = NaN(length(RandomEffects)-1, nyVars); 
ivec         = 1:nyVars;
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
allSeeds = randi(999999, length(nObservations), 1);

%% Summary of all ground truth and settings
settings.nObservations          = nObservations;
settings.nXvars                 = nXvars;
settings.nyVars                 = nyVars;
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

save(fullfile(outDir, 'Settings.mat'), 'settings');

%% Start experiment
for obs = 1:length(nObservations)
    
    % Current number of observations
    currObs = nObservations(obs);
    
    % Initialize a few variables
    eid              = ones(currObs, 1);
    agevec           = zeros(currObs, 1);
    y_RFX            = nan(currObs, nyVars);
    
    % Set seed for this nObservation value
    rng(allSeeds(obs), 'twister');

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
    
    %% Additionally generate a pihatmat
    % pihatmat is a square matrix with size matching the number of 
    % individuals (not the number of observations)
    [~, IA]         = unique(iid, 'stable');
    pihatmat        = eye(length(IA));
    [~, ~, IC_fam]  = unique(fid(IA), 'stable'); 
    nfam            = max(IC_fam);
    for fi          = 1:nfam
        idxvec      = rowvec(find(IC_fam==fi)); % Identify all subjects for a given family
        t           = triu(0.5 + 0.05 * randn(length(idxvec)), 1);
        pihatmat(idxvec, idxvec) = pihatmat(idxvec, idxvec) + t + t';
    end

    % Using FASE model, generate y_RFX
    % Parse family structure
    clusterinfo = FEMA_parse_family(iid, eid, fid, agevec, pihatmat, ...
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
    
    % Output directory
    tmpOutDir = fullfile(outDir, ['nObs-', num2str(currObs, '%04d')]);
    if not(exist(tmpOutDir, 'dir'))
        mkdir(tmpOutDir);
    end
    
    %% FE models
    % Initialize timer for FEMA
    FEMA_FE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_FE, beta_se_FE, zmat_FE, logpmat_FE, sig2tvec_FE, sig2mat_FE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),          ...
              allBins, [], 'RandomEffects', {'F', 'E'});
            
    % End timer for FEMA
    FEMA_FE_elapsed = toc(FEMA_FE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_FE.mat'), 'beta_hat_FE',    'beta_se_FE', 'zmat_FE', 'logpmat_FE', 'sig2tvec_FE', 'sig2mat_FE', ...
                                                     'FEMA_FE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                         
    % Now loop over all phenotypes and estimate model using fitlmematrix
    mdls_FE      = cell(nyVars, 1);
    lme_FE_Init  = tic;
    for yVar = 1:nyVars
        mdls_FE{yVar} = fitlmematrix(X, ymat(:,yVar), {ones(length(fid),1)}, {fid'});
    end
    lme_FE_Elapsed = toc(lme_FE_Init);
    
    % Extract out everything relevant
    results_FE = table('Size', [nyVars, 20], 'VariableTypes', repmat({'double'}, 20, 1), 'VariableNames', header);
    for yVar = 1:nyVars
        [a, b, c]                  = covarianceParameters(mdls_FE{yVar});
        results_FE.YvarNum(yVar)   = yVar;
        results_FE.Beta_X1(yVar)   = mdls_FE{yVar}.Coefficients.Estimate(1);
        results_FE.Beta_X2(yVar)   = mdls_FE{yVar}.Coefficients.Estimate(2);
        results_FE.Beta_X3(yVar)   = mdls_FE{yVar}.Coefficients.Estimate(3);
        results_FE.Beta_X4(yVar)   = mdls_FE{yVar}.Coefficients.Estimate(4);
        results_FE.Beta_X5(yVar)   = mdls_FE{yVar}.Coefficients.Estimate(5);
        results_FE.SE_X1(yVar)     = mdls_FE{yVar}.Coefficients.SE(1);
        results_FE.SE_X2(yVar)     = mdls_FE{yVar}.Coefficients.SE(2);
        results_FE.SE_X3(yVar)     = mdls_FE{yVar}.Coefficients.SE(3);
        results_FE.SE_X4(yVar)     = mdls_FE{yVar}.Coefficients.SE(4);
        results_FE.SE_X5(yVar)     = mdls_FE{yVar}.Coefficients.SE(5);
        results_FE.V_F(yVar)       = c{1}.Estimate;
        results_FE.V_E(yVar)       = c{2}.Estimate;
        results_FE.Low_VF(yVar)    = c{1}.Lower;
        results_FE.Low_VE(yVar)    = c{2}.Lower;
        results_FE.Upp_VF(yVar)    = c{1}.Upper;
        results_FE.Upp_VE(yVar)    = c{2}.Upper;
    end
    
%     % Prepare table of values first
%     % Ignore the first X column as intercept
%     tblData      = [fid', iid', num2cell([X(:,2:end), ymat])];
%     tblNames     = [{'FID', 'IID', 'X2', 'X3', 'X4', 'X5'}, strcat({'y'}, num2str((1:nyVars)', '%02d'))'];
%     tblData      = cell2table(tblData, 'VariableNames', tblNames);
%     mdls_FE      = cell(nyVars, 1);
%     lme_FE_Init  = tic;
%     for yVar = 1:nyVars
%         mdls_FE{yVar} = fitlme(tblData, ['y', num2str(yVar, '%02d'), ' ~ 1 + X2 + X3 + X4 + X5 + (1|FID)'], 'FitMethod', 'ML');
%     end
%     lme_FE_Elapsed = toc(lme_FE_Init);
    
    % Save everything
    save(fullfile(tmpOutDir, 'Results_fitlme_FE.mat'), 'results_FE', 'lme_FE_Elapsed');
    
    % Update user
    disp(['Finished: FE: Sample size ', num2str(currObs, '%04d')]);
    
    %% SE models
    % Initialize timer for FEMA
    FEMA_SE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_SE, beta_se_SE, zmat_SE, logpmat_SE, sig2tvec_SE, sig2mat_SE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),          ...
              allBins, [], 'RandomEffects', {'S', 'E'});
            
    % End timer for FEMA
    FEMA_SE_elapsed = toc(FEMA_SE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_SE.mat'), 'beta_hat_SE',    'beta_se_SE', 'zmat_SE', 'logpmat_SE', 'sig2tvec_SE', 'sig2mat_SE', ...
                                                    'FEMA_SE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                                
    % Now loop over all phenotypes and estimate model using fitlmematrix
    mdls_SE      = cell(nyVars, 1);
    lme_SE_Init  = tic;
    for yVar = 1:nyVars
        mdls_SE{yVar} = fitlmematrix(X, ymat(:,yVar), {ones(length(iid),1)}, {iid'});
    end
    lme_SE_Elapsed = toc(lme_SE_Init);
    
    % Extract out everything relevant
    results_SE = table('Size', [nyVars, 20], 'VariableTypes', repmat({'double'}, 20, 1), 'VariableNames', header);
    for yVar = 1:nyVars
        [a, b, c]                   = covarianceParameters(mdls_SE{yVar});
        results_SE.YvarNum(yVar)   = yVar;
        results_SE.Beta_X1(yVar)   = mdls_SE{yVar}.Coefficients.Estimate(1);
        results_SE.Beta_X2(yVar)   = mdls_SE{yVar}.Coefficients.Estimate(2);
        results_SE.Beta_X3(yVar)   = mdls_SE{yVar}.Coefficients.Estimate(3);
        results_SE.Beta_X4(yVar)   = mdls_SE{yVar}.Coefficients.Estimate(4);
        results_SE.Beta_X5(yVar)   = mdls_SE{yVar}.Coefficients.Estimate(5);
        results_SE.SE_X1(yVar)     = mdls_SE{yVar}.Coefficients.SE(1);
        results_SE.SE_X2(yVar)     = mdls_SE{yVar}.Coefficients.SE(2);
        results_SE.SE_X3(yVar)     = mdls_SE{yVar}.Coefficients.SE(3);
        results_SE.SE_X4(yVar)     = mdls_SE{yVar}.Coefficients.SE(4);
        results_SE.SE_X5(yVar)     = mdls_SE{yVar}.Coefficients.SE(5);
        results_SE.V_S(yVar)       = c{1}.Estimate;
        results_SE.V_E(yVar)       = c{2}.Estimate;
        results_SE.Low_VS(yVar)    = c{1}.Lower;
        results_SE.Low_VE(yVar)    = c{2}.Lower;
        results_SE.Upp_VS(yVar)    = c{1}.Upper;
        results_SE.Upp_VE(yVar)    = c{2}.Upper;
    end

%     % Now loop over all phenotypes and estimate model using fitlme
%     % Prepare table of values first
%     % Ignore the first X column as intercept
%     tblData      = [fid', iid', num2cell([X(:,2:end), ymat])];
%     tblNames     = [{'FID', 'IID', 'X2', 'X3', 'X4', 'X5'}, strcat({'y'}, num2str((1:nyVars)', '%02d'))'];
%     tblData      = cell2table(tblData, 'VariableNames', tblNames);
%     mdls_SE      = cell(nyVars, 1);
%     lme_SE_Init  = tic;
%     for yVar = 1:nyVars
%         mdls_SE{yVar} = fitlme(tblData, ['y', num2str(yVar, '%02d'), ' ~ 1 + X2 + X3 + X4 + X5 + (1|IID)'], 'FitMethod', 'ML');
%     end
%     lme_SE_Elapsed = toc(lme_FE_Init);
    
    % Save everything
    save(fullfile(tmpOutDir, 'Results_fitlme_SE.mat'), 'results_SE', 'lme_SE_Elapsed');
    
    % Update user
    disp(['Finished: SE: Sample size ', num2str(currObs, '%04d')]);
    
    %% FSE models
    % Initialize timer for FEMA
    FEMA_FSE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_FSE, beta_se_FSE, zmat_FSE, logpmat_FSE, sig2tvec_FSE, sig2mat_FSE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),          ...
              allBins, [], 'RandomEffects', {'F', 'S', 'E'});
            
    % End timer for FEMA
    FEMA_FSE_elapsed = toc(FEMA_FSE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_FSE.mat'), 'beta_hat_FSE',    'beta_se_FSE', 'zmat_FSE', 'logpmat_FSE', 'sig2tvec_FSE', 'sig2mat_FSE', ...
                                                      'FEMA_FSE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                         
    % Now loop over all phenotypes and estimate model using fitlmematrix
    mdls_FSE      = cell(nyVars, 1);
    lme_FSE_Init  = tic;
    for yVar = 1:nyVars
        mdls_FSE{yVar} = fitlmematrix(X, ymat(:,yVar), {ones(length(fid),1), ones(length(iid),1)}, {fid', iid'});
    end
    lme_FSE_Elapsed = toc(lme_FSE_Init);
    
    % Extract out everything relevant
    results_FSE = table('Size', [nyVars, 20], 'VariableTypes', repmat({'double'}, 20, 1), 'VariableNames', header);
    for yVar = 1:nyVars
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
    %     % Now loop over all phenotypes and estimate model using fitlme
%     % Prepare table of values first
%     % Ignore the first X column as intercept
%     tblData      = [fid', iid', num2cell([X(:,2:end), ymat])];
%     tblNames     = [{'FID', 'IID', 'X2', 'X3', 'X4', 'X5'}, strcat({'y'}, num2str((1:nyVars)', '%02d'))'];
%     tblData      = cell2table(tblData, 'VariableNames', tblNames);
%     mdls_FSE     = cell(nyVars, 1);
%     lme_FSE_Init = tic;
%     for yVar = 1:nyVars
%         mdls_FSE{yVar} = fitlme(tblData, ['y', num2str(yVar, '%02d'), ' ~ 1 + X2 + X3 + X4 + X5 + (1|FID) + (1|IID)'], 'FitMethod', 'ML');
%     end
%     lme_FSE_Elapsed = toc(lme_FSE_Init);
    
    % Save everything
    save(fullfile(tmpOutDir, 'Results_fitlme_FSE.mat'), 'results_FSE', 'lme_FSE_Elapsed');
    
    % Update user
    disp(['Finished: FSE: Sample size ', num2str(currObs, '%04d')]);
    
    %% Additionally do models that include A term but not for MATLAB
    %% AE models
    % Initialize timer for FEMA
    FEMA_AE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_AE, beta_se_AE, zmat_AE, logpmat_AE, sig2tvec_AE, sig2mat_AE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),          ...
              allBins, pihatmat, 'RandomEffects', {'A', 'E'});
            
    % End timer for FEMA
    FEMA_AE_elapsed = toc(FEMA_AE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_AE.mat'), 'beta_hat_AE',    'beta_se_AE', 'zmat_AE', 'logpmat_AE', 'sig2tvec_AE', 'sig2mat_AE', ...
                                                     'FEMA_AE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                         
    % Update user
    disp(['Finished: AE: Sample size ', num2str(currObs, '%04d')]);
    
    %% FAE models
    % Initialize timer for FEMA
    FEMA_FAE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_FAE, beta_se_FAE, zmat_FAE, logpmat_FAE, sig2tvec_FAE, sig2mat_FAE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),          ...
              allBins, pihatmat, 'RandomEffects', {'F', 'A', 'E'});
            
    % End timer for FEMA
    FEMA_FAE_elapsed = toc(FEMA_FAE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_FAE.mat'), 'beta_hat_FAE',    'beta_se_FAE', 'zmat_FAE', 'logpmat_FAE', 'sig2tvec_FAE', 'sig2mat_FAE', ...
                                                     'FEMA_FAE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                         
    % Update user
    disp(['Finished: FAE: Sample size ', num2str(currObs, '%04d')]);
    
    %% SAE models
    % Initialize timer for FEMA
    FEMA_SAE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_SAE, beta_se_SAE, zmat_SAE, logpmat_SAE, sig2tvec_SAE, sig2mat_SAE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),          ...
              allBins, pihatmat, 'RandomEffects', {'S', 'A', 'E'});
            
    % End timer for FEMA
    FEMA_SAE_elapsed = toc(FEMA_SAE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_SAE.mat'), 'beta_hat_SAE',    'beta_se_SAE', 'zmat_SAE', 'logpmat_SAE', 'sig2tvec_SAE', 'sig2mat_SAE', ...
                                                     'FEMA_SAE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                         
    % Update user
    disp(['Finished: SAE: Sample size ', num2str(currObs, '%04d')]);
    
    %% FASE models
    % Initialize timer for FEMA
    FEMA_FASE_timeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat_FASE, beta_se_FASE, zmat_FASE, logpmat_FASE, sig2tvec_FASE, sig2mat_FASE] = ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),          ...
              allBins, pihatmat, 'RandomEffects', {'F', 'A', 'S', 'E'});
            
    % End timer for FEMA
    FEMA_FASE_elapsed = toc(FEMA_FASE_timeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA_FASE.mat'), 'beta_hat_FASE',    'beta_se_FASE', 'zmat_FASE', 'logpmat_FASE', 'sig2tvec_FASE', 'sig2mat_FASE', ...
                                                     'FEMA_FASE_elapsed', 'ymat', 'X', 'iid', 'fid');
                                         
    % Update user
    disp(['Finished: FASE: Sample size ', num2str(currObs, '%04d')]);
    
    %% Initialize parallel pool
    % Start a pool
    local            = parcluster('local');
    local.NumThreads = 2;
    pool             = local.parpool(20, 'IdleTimeout', 240);
    
    %% Other initialization
    mdls_parFE       = cell(nyVars, 1);
    mdls_parSE       = cell(nyVars, 1);
    mdls_parFSE      = cell(nyVars, 1);
    
    %% Redo processing with parallel computing - FE
    lme_parFE_Init   = tic;
    parfor yVar = 1:nyVars
        mdls_parFE{yVar} = fitlmematrix(X, ymat(:,yVar), {ones(length(fid),1)}, {fid'});
    end
    lme_parFE_Elapsed = toc(lme_parFE_Init);
    
    %% SE model
    lme_parSE_Init   = tic;
    parfor yVar = 1:nyVars
        mdls_parSE{yVar} = fitlmematrix(X, ymat(:,yVar), {ones(length(iid),1)}, {iid'});
    end
    lme_parSE_Elapsed = toc(lme_parSE_Init);

    %% FSE model
    lme_parFSE_Init   = tic;
    parfor yVar = 1:nyVars
        mdls_parFSE{yVar} = fitlmematrix(X, ymat(:,yVar), {ones(length(fid),1), ones(length(iid),1)}, {fid', iid'});
    end
    lme_parFSE_Elapsed = toc(lme_parFSE_Init);
    
    % Delete the pool
    delete(pool);
    
    % Save timing
    save(fullfile(tmpOutDir, 'Results_fitlme_parTimes.mat'), 'lme_parFE_Elapsed', 'lme_parSE_Elapsed', 'lme_parFSE_Elapsed');
end