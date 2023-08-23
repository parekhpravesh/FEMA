%% Simulation 002: comparison with fitlme
%% Include relevant paths
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/cmig_tools_internal-beta/cmig_tools_utils/matlab');
addpath('/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/cmig_tools_internal-beta/FEMA');

%% Prepare output directory
outDir = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/Simulation002_Compare_fitlme_perms_nn';
if not(exist(outDir, 'dir'))
    mkdir(outDir);
end

%% Setup parameters 
nObservations    = linspace(2000, 20000, 10);
nXvars           = 5;
nyVars           = 50;
nFamMembers      = 5;
nRepObservations = 5;
RandomEffects    = {'F', 'S', 'E'}; 
epsmin           = 0.2; 
sigmaLow         = 0.2;
allBins          = 20;
niter            = 1;
numPerms         = 100;

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
    tmpOutDir = fullfile(outDir, ['nObs-', num2str(currObs, '%04d')]);
    if not(exist(tmpOutDir, 'dir'))
        mkdir(tmpOutDir);
    end
        
    % Initialize timer for FEMA
    FEMAtimeInit = tic;

    % Eatimate model parameters using FEMA
    [beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat,       ...
     Hessmat, logLikvec, beta_hat_perm, beta_se_perm,           ...
     zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm] =  ...
     FEMA_fit(X, iid, eid, fid, agevec, ymat, niter, ones(1,nXvars),  ...
              allBins, [], 'RandomEffects', RandomEffects, 'nperms', numPerms, 'PermType', 'wildbootstrap-nn');
            
    % End timer for FEMA
    FEMAelapsed = toc(FEMAtimeInit);
        
    % Save everything
    save(fullfile(tmpOutDir, 'Results_FEMA.mat'), 'beta_hat',    'beta_se', 'zmat', 'logpmat', 'sig2tvec', 'sig2mat',     ...
                                                  'FEMAelapsed', 'Hessmat', 'logLikvec', 'beta_hat_perm', 'beta_se_perm', ...
                                                  'zmat_perm',   'sig2tvec_perm', 'sig2mat_perm', 'logLikvec_perm', 'ymat', 'X', 'iid', 'fid');
                                         
    % Now loop over all phenotypes and estimate model using fitlme
    % Prepare table of values first
    % Ignore the first X column as intercept
    tblData  = [fid', iid', num2cell([X(:,2:end), ymat])];
    tblNames = [{'FID', 'IID', 'X2', 'X3', 'X4', 'X5'}, strcat({'y'}, num2str((1:nyVars)', '%02d'))'];
    tblData  = cell2table(tblData, 'VariableNames', tblNames);
    mdls     = cell(nyVars, 1);
    lmeInit  = tic;
    for yVar = 1:nyVars
        mdls{yVar} = fitlme(tblData, ['y', num2str(yVar, '%02d'), ' ~ 1 + X2 + X3 + X4 + X5 + (1|FID) + (1|IID)'], 'FitMethod', 'ML');
    end
    lmeElapsed = toc(lmeInit);
    
    % Save everything
    save(fullfile(tmpOutDir, 'Results_fitlme.mat'), 'mdls', 'lmeElapsed', 'tblData');
    
    % Update user
    disp(['Finished: Sample size ', num2str(currObs, '%04d')]);
end