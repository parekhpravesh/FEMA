%% Simulation 005: comparison with fitlme as a function of number of observations
%% Include relevant paths
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/codebase/cmig_tools-2.3.0/cmig_tools_utils/matlab');
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/codebase/cmig_tools-2.3.0/FEMA');

%% Prepare output directory
outDir = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/Simulation005_Compare_fitlme_perms_nn_nObs';
if not(exist(outDir, 'dir'))
    mkdir(outDir);
end

%% Setup parameters 
nObservations    = [50, 100:100:1000, 2000:2000:10000];
nXvars           = 5;
nyVars           = 500;
nFamMembers      = 5;
nRepObservations = 5;
RandomEffects    = {'F', 'S', 'E'}; 
epsmin           = 0.2; 
sigmaLow         = 0.2;
allBins          = 20;
niter            = 1;
numPerms         = 100;
numRepeats       = 5;

%% Generate seeds
rng(20230809, 'twister');
allSeeds = unique(randi(999999, length(nObservations) * numRepeats * 2, 1), 'stable');
count    = 1;

%% Generate ground truth effect sizes
rng(allSeeds(count), 'twister');

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

% Summary of all ground truth and settings
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

count = count + 1;

%% Start experiment
for obs = 1:length(nObservations)
    % Current number of observations
    currObs = nObservations(obs);
    
    % Initialize a few variables
    eid              = ones(currObs, 1);
    agevec           = zeros(currObs, 1);
    y_RFX            = nan(currObs, nyVars);
    
    for repeats = 1:numRepeats
    
        % Set seed for this nObservation value
        rng(allSeeds(count), 'twister');
        
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
                            sprintf('V_%s',RandomEffects{ri}))),nyVars)'; %#ok<GFLD>
            end
            y_RFX(clusterinfo{fi}.jvec_fam,:) = sqrt(sig2tvec_true) .* tmp;
        end

        % Put together as a single y variable
        ymat = y_FFX + y_RFX;

        % Output directory
        tmpOutDir = fullfile(outDir, [num2str(obs, '%04d'), '-nObs-', num2str(currObs, '%04d')], [num2str(obs, '%04d'), '-nObs-', num2str(currObs, '%04d'), '-Repeat-', num2str(repeats, '%02d')]);
        if not(exist(tmpOutDir, 'dir'))
            mkdir(tmpOutDir);
        end
        
        % Initialize timer for FEMA
        FEMAtimeInit = tic;

        % Eatimate model parameters using FEMA
        [beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat,               ...
         Hessmat, logLikvec, beta_hat_perm, beta_se_perm,                   ...
         zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm,            ...
         binvec_save, nvec_bins, tvec_bins, FamilyStruct, reusableVars] =   ...
         FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),    ...
                  allBins, [], 'RandomEffects', RandomEffects,              ...
                  'nperms', numPerms, 'PermType', 'wildbootstrap-nn', 'returnReusable', true);
            
        % End timer for FEMA
        FEMAelapsed = toc(FEMAtimeInit);

        % Save everything
        save(fullfile(tmpOutDir, 'Results_FEMA.mat'), 'beta_hat',    'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat',       ...
                                                      'FEMAelapsed', 'Hessmat', 'logLikvec', 'beta_hat_perm', 'beta_se_perm',   ...
                                                      'zmat_perm',   'sig2tvec_perm', 'sig2mat_perm', 'logLikvec_perm', 'ymat', ...
                                                      'X', 'iid', 'fid', 'binvec_save', 'FamilyStruct', 'reusableVars');
                                         
        % Now loop over all phenotypes and estimate model using fitlme
        % Call using fitlmematrix so no time is spent on making a table type
        % variable
        % Use parallel computing to speed up things
        % Start a pool
        local            = parcluster('local');
        local.NumThreads = 2;
        pool             = local.parpool(10, 'IdleTimeout', 240);

        mdls             = cell(nyVars, 1);
        Zvariable        = {ones(length(fid),1), ones(length(iid),1)};
        Gvariable        = {fid', iid'};
        lmeInit          = tic;
        parfor yVar      = 1:nyVars
            mdls{yVar}   = fitlmematrix(X, ymat(:,yVar), Zvariable, Gvariable);
        end
        lmeElapsed       = toc(lmeInit);

        % Delete the pool
        delete(pool);
        
        % Extract out LME parameters 
        [beta_hat_fitlme, beta_se_fitlme]                   = deal(zeros(settings.nXvars, settings.nyVars));
        [sig2mat_fitlme,  sig2Low_fitlme, sig2Upp_fitlme]   = deal(zeros(3, settings.nyVars));
        
        for nY = 1:settings.nyVars
            [psi, mse, stats]                      = covarianceParameters(mdls{nY,1});
            beta_hat_fitlme(1:settings.nXvars, nY) = mdls{nY,1}.Coefficients.Estimate;
            beta_se_fitlme(1:settings.nXvars,  nY) = mdls{nY,1}.Coefficients.SE;
            sig2mat_fitlme(1:3,                nY) = [stats{1}.Estimate; stats{2}.Estimate; stats{3}.Estimate];
            sig2Low_fitlme(1:3,                nY) = [stats{1}.Lower; stats{2}.Lower; stats{3}.Lower];
            sig2Upp_fitlme(1:3,                nY) = [stats{1}.Upper; stats{2}.Upper; stats{3}.Upper];
        end
        
        % Remember that these are in standard deviation units - to be squared
        % before comparing with FEMA output
        variables_lme.beta_hat_fitlme = beta_hat_fitlme;
        variables_lme.beta_se_fitlme  = beta_se_fitlme;
        variables_lme.sig2mat_fitlme  = sig2mat_fitlme;
        variables_lme.sig2Low_fitlme  = sig2Low_fitlme;
        variables_lme.sig2Upp_fitlme  = sig2Upp_fitlme;
        variables_lme.elapsed         = lmeElapsed;
    
        % Save
        % No longer saving each model - takes a huge amount of resources
        % save(fullfile(tmpOutDir, 'Results_fitlmematrix.mat'), 'mdls', 'lmeElapsed', '-v7.3');
        save(fullfile(tmpOutDir, 'Results_fitlmematrix_summary.mat'), 'variables_lme');
    
        % Update user
        disp(['Finished: Sample size ', num2str(currObs, '%04d'), ', Repeat ', num2str(repeats, '%02d')]);
        
        % Update count
        count = count + 1;
    end
end