%% Analyze sMRI DK cortical thickness measures vertex wise
%% Setup
% Add FEMA code to path
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/cmig_tools-2.3.0/cmig_tools_utils/matlab');
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/cmig_tools-2.3.0/showSurf');
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/cmig_tools-2.3.0/FEMA');

% Specify release number
dataRelease = '4.0';

% Various paths
dirABCD         = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/';
dirOut          = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/CorticalThickness_vertexWise_FSE';
dirSupport      = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered/';
dirDesign       = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/';
dirTabulated    = fullfile(dirABCD, dataRelease, 'tabulated', 'released'); 
dirImaging      = fullfile(dirABCD, dataRelease, 'imaging_concat', 'vertexwise', 'smri');
fnameGRM        = fullfile(dirABCD, dataRelease, 'genomics', ['ABCD_rel', dataRelease, '_grm.mat']);
fnameDesign     = 'GlobalDesignMatrix.csv';
fnameImaging    = 'thickness_ic5_sm256';

% Prepare output directory
if ~exist(dirOut, 'dir')
    mkdir(dirOut);
end

%% Read various files
% First, read the QC variables
dataQC = readtable(fullfile(dirABCD, dataRelease, 'tabulated', 'img', 'abcd_imgincl01.csv'));

% Next, get some additional info about MRI acquisition
dataInfo = readtable(fullfile(dirABCD, dataRelease, 'tabulated', 'img', 'abcd_mri01.csv'));

% Read various support files - only need allPCs and imgVars
imgVars       = readtable(fullfile(dirSupport, 'ABCD_rel4.0_covars_img_base_2yr.txt'));
allPCs        = readtable(fullfile(dirSupport, 'ABCD_rel4.0_pcs_base_2yr.txt'));

% Unique event names - useful later
allUqEvents = unique(dataQC.eventname);

%% Next, read the design matrix
designMatrix    = readtable(fullfile(dirDesign, fnameDesign));

%% Get vertex-wise data
[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat] = FEMA_process_data(fnameImaging, dirTabulated, dirImaging, 'vertex', 'pihat_file', fnameGRM);

%% Within the design matrix, remove subjects who did not pass QC
toDelete = unique(designMatrix.src_subject_id(designMatrix.imgincl_t1w_include == 0));
designMatrix(ismember(designMatrix.src_subject_id, toDelete), :) = [];

%% Add a column for alignment
ymat_toAlign = strcat(iid_concat, eid_concat);

% Ensure that data lines up with design matrix
[a, ~]       = ismember(ymat_toAlign, designMatrix.toAlign);
ymat         = ymat(a, :);
iid_concat   = iid_concat(a, :);
eid_concat   = eid_concat(a, :);
ymat_toAlign = ymat_toAlign(a, :);

[a, b]       = ismember(designMatrix.toAlign, ymat_toAlign);
ymat         = ymat(b, :);
iid_concat   = iid_concat(b, :);
eid_concat   = eid_concat(b, :);
ymat_toAlign = ymat_toAlign(b, :);

%% Attempt to create covariates
% To use: age at recruitment, age delta, sex, scanner, household income, 
% highest parental education, and the first 20 genetic PCs

% Location for PCs
locGenesisPCs = ~cellfun(@isempty, regexpi(designMatrix.Properties.VariableNames, '^genesis'));

intercept       = ones(height(designMatrix), 1);
geneticPCs      = designMatrix{:, locGenesisPCs};

% Compute ageRecruitment - same as baseline age
locBaseline     = strcmpi(designMatrix.eventname, 'baseline_year_1_arm_1');
ageRecruitment  = [designMatrix.interview_age(locBaseline); designMatrix.interview_age(locBaseline)];

% Compute ageDelta
ageDelta                 = zeros(height(designMatrix), 1);
ageDelta(~locBaseline,1) = designMatrix.interview_age(~locBaseline) - designMatrix.interview_age(locBaseline);

% Sex
sex                  = categorical(designMatrix.sex);
namesSex             = categories(sex);
varSex               = dummyvar(sex);

% MRI device number
deviceInfo           = categorical(designMatrix.mri_info_deviceserialnumber);
namesDeviceInfo      = categories(deviceInfo);
varDeviceInfo        = dummyvar(deviceInfo);

% MRI software number
softwareInfo         = categorical(designMatrix.mri_info_softwareversion);
namesSoftwareInfo    = categories(softwareInfo);
varSoftwareInfo      = dummyvar(softwareInfo);

% Household income level
householdIncome      = categorical(designMatrix.household_income);
namesHouseholdIncome = categories(householdIncome);
varHouseholdIncome   = dummyvar(householdIncome);

% Parental educational level
parentalEdu             = categorical(designMatrix.high_educ);
namesParentalEducation  = categories(parentalEdu);
varParentalEducation    = dummyvar(parentalEdu);

% Put covariates together and put covariate names together
% Only keep the first level of sex
% Keep n-2 level of device info
% Keep n-2 level of software info
% Keep n-1 level of household income
% Keep n-1 level of parental education level
% Keeping n-1 levels of device info and software info leads to rank
% deficient matrix - why? Likely because there are too few observations
covariates = [intercept, ageRecruitment, ageDelta, varSex(:,1),         ...
             varDeviceInfo(:, 1:end-2), varSoftwareInfo(:, 1:end-2),    ...
             varHouseholdIncome(:, 1:end-1), varParentalEducation(:, 1:end-1), geneticPCs];
         
%% Ensure no rank deficiency
if rank(covariates) ~= size(covariates,2)
    error('Rank deficient covariates matrix; check covariates');
end
         
% Names of covariates
names_genesis = cellfun(@(x) strrep(x, ' ', ''), strcat({'genensis_PC'}, num2str((1:20)')), 'UniformOutput', false);
covarNames    = [{'Intercept', 'ageBaseline', 'ageDelta'}, namesSex(1)', namesDeviceInfo(1:end-2)', ...
                 namesSoftwareInfo(1:end-2)', namesHouseholdIncome(1:end-1)', namesParentalEducation(1:end-1)', names_genesis'];
             
%% Create a new version of the design matrix with covariates
% First four columns should be: src_subject_id, eventname, rel_family_id, age
toIntersect = cell2table([designMatrix.src_subject_id, designMatrix.eventname, num2cell([designMatrix.rel_family_id, designMatrix.interview_age, covariates])], 'VariableNames', ...
                         [{'src_subject_id', 'eventname', 'rel_family_id', 'age'}, covarNames]);

% Save as T1w design matrix
writetable(toIntersect, fullfile(dirOut, 'DesignMatrix_T1w_CorticalThickness.csv'));

%% Intersect the design matrix and the vertex-wise data
[X, iid, eid, fid, agevec, ymatUse, contrasts, colnames_model, pihatmat, PregID, HomeID] = FEMA_intersect_design(fullfile(dirOut, 'DesignMatrix_T1w_CorticalThickness.csv'), ...
                                                                                                              ymat, iid_concat, eid_concat, 'pihat', pihat);

%% Standardize ymatUse
ymatUse_std = (ymatUse - mean(ymatUse))./std(ymatUse);

%% Prepare random effects
SubjectEffect = iid;
FamilyEffect  = cellfun(@(x) strrep(x, ' ', ''), strcat({'F'}, num2str(fid)), 'UniformOutput', false);

tempVar = [designMatrix.interview_age(locBaseline), designMatrix.interview_age(~locBaseline)];

%% Additional variables for FEMA_fit
nperms          = 0;
RandomEffects   = {'F', 'S', 'E'};
niter           = 1;
nbins           = 20;

%% Call FEMA_fit
initFEMA = tic;
[beta_hat,      beta_se,        zmat,        logpmat,       ...
 sig2tvec,      sig2mat,        Hessmat,     logLikvec,     ...
 beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm, ...
 sig2mat_perm,  logLikvec_perm, binvec_save] =              ...
 FEMA_fit(X, SubjectEffect, eid, FamilyEffect, agevec, ymatUse_std, niter,    ...
          contrasts, nbins, [], 'RandomEffects', RandomEffects, 'nperms', nperms);
elapsedFEMA = toc(initFEMA);

%% Save everything
if not(exist(dirOut, 'dir'))
    mkdir(dirOut);
end
save(fullfile(dirOut, 'Results.mat'), '-v7.3');

%% Additional saving borrwed from FEMA_wrapper
if sum(~mask)>0
    
    z_tmp        = zeros(size(zmat,1),size(mask,2));
    p_tmp        = zeros(size(logpmat,1),size(mask,2));
    beta_tmp     = zeros(size(beta_hat,1),size(mask,2));
    betase_tmp   = zeros(size(beta_se,1),size(mask,2));
    sig2mat_tmp  = zeros(size(sig2mat,1),size(mask,2));
    sig2tvec_tmp = zeros(size(sig2tvec,1),size(mask,2));
    
    z_tmp(:,ivec_mask)       = zmat;
    p_tmp(:,ivec_mask)       = logpmat;
    beta_tmp(:,ivec_mask)    = beta_hat;
    betase_tmp(:,ivec_mask)  = beta_se;
    sig2mat_tmp(:,ivec_mask) = sig2mat;
    sig2tvec_tmp(:,ivec_mask)= sig2tvec;
    
    zmat     = z_tmp;
    logpmat  = p_tmp;
    beta_hat = beta_tmp;
    beta_se  = betase_tmp;
    sig2mat  = sig2mat_tmp;
    sig2tvec = sig2tvec_tmp;
    
    if nperms > 0
        zperm_tmp           = zeros(size(zmat_perm,1),size(mask,2),size(zmat_perm,3));
        betaperm_tmp        = zeros(size(beta_hat_perm,1),size(mask,2),size(beta_hat_perm,3));
        betaseperm_tmp      = zeros(size(beta_se_perm,1),size(mask,2),size(beta_se_perm,3));
        sig2matperm_tmp     = zeros(size(sig2mat_perm,1),size(mask,2),size(sig2mat_perm,3));
        sig2tvecperm_tmp    = zeros(size(sig2tvec_perm,1),size(mask,2),size(sig2tvec_perm,3));
        
        zperm_tmp(:,ivec_mask,:)        = zmat_perm;
        betaperm_tmp(:,ivec_mask,:)     = beta_hat_perm;
        betaseperm_tmp(:,ivec_mask,:)   = beta_se_perm;
        sig2matperm_tmp(:,ivec_mask,:)  = sig2mat_perm;
        sig2tvecperm_tmp(:,ivec_mask,:) = sig2tvec_perm;
        
        zmat_perm       = zperm_tmp;
        beta_hat_perm   = betaperm_tmp;
        beta_se_perm    = betaseperm_tmp;
        sig2mat_perm    = sig2matperm_tmp;
        sig2tvec_perm   = sig2tvecperm_tmp; 
    end
end

% Save output
base_variables_to_save = {'X','iid','eid','colnames_model','contrasts','zmat','logpmat','beta_hat','beta_se','sig2mat','sig2tvec','mask', 'binvec_save'};
save(fullfile(dirOut, 'Results_masked.mat'), base_variables_to_save{:}, '-v7.3');

%% Generate lookup tables
fsdir                       = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/from_RH4-x86_64-R711';
cd(fsdir);
lookup_age_uncorrected      = FEMA_lookupVertices(beta_hat(2,:), zmat(2,:), beta_se(2,:), [], [],           fsdir);
lookup_age_Bonferroni       = FEMA_lookupVertices(beta_hat(2,:), zmat(2,:), beta_se(2,:), [], 0.05/18742,   fsdir);
lookup_ageDelta_uncorrected = FEMA_lookupVertices(beta_hat(3,:), zmat(3,:), beta_se(3,:), [], [],           fsdir);
lookup_ageDelta_Bonferroni  = FEMA_lookupVertices(beta_hat(3,:), zmat(3,:), beta_se(3,:), [], 0.05/18742,   fsdir);

writetable(lookup_age_uncorrected{1},       fullfile(dirOut, 'Lookup_age_uncorrected.xlsx'));
writetable(lookup_age_Bonferroni{1},        fullfile(dirOut, 'Lookup_age_Bonferroni.xlsx'));
writetable(lookup_ageDelta_uncorrected{1},  fullfile(dirOut, 'Lookup_ageDelta_uncorrected.xlsx'));
writetable(lookup_ageDelta_Bonferroni{1},   fullfile(dirOut, 'Lookup_ageDelta_Bonferroni.xlsx'));