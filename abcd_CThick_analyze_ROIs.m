%% Analyze sMRI DK cortical thickness measures across ROIs using FSE model (rate of change) and compare with MATLAB
%% Setup
% Add FEMA code to path
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/cmig_tools-2.3.0/cmig_tools_utils/matlab');
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/cmig_tools-2.3.0/FEMA');

% Specify release number
dataRelease = '4.0';

% Various paths
dirABCD         = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/';
dirOut          = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/CorticalThickness_DK_FSE_nn';
dirSupport      = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered/';
dirDesign       = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone';
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

%% Read some cortical thickness data across ROIs
dataCThick = readtable(fullfile(dirABCD, dataRelease, 'tabulated', 'img', 'abcd_smrip102.csv'));

%% Within the design matrix, remove subjects who did not pass QC
toDelete = unique(designMatrix.src_subject_id(designMatrix.imgincl_t1w_include == 0));
designMatrix(ismember(designMatrix.src_subject_id, toDelete), :) = [];

%% Add a column for table alignment
dataCThick.toAlign = strcat(dataCThick.src_subject_id, dataCThick.eventname);

%% Ensure that cortical thickness data lines up with design matrix
[a, ~]      = ismember(dataCThick.toAlign, designMatrix.toAlign);
dataCThick  = dataCThick(a, :);
[a, b]      = ismember(designMatrix.toAlign, dataCThick.toAlign);
dataCThick  = dataCThick(b, :);

%% Identify ROIs of interest
varsInterest = ~cellfun(@isempty, regexpi(dataCThick.Properties.VariableNames, '^smri_thick_cdk'));
yVariables = dataCThick(:, varsInterest);
yVariables = yVariables(:, cellfun(@isempty, regexpi(yVariables.Properties.VariableNames, '1$')));
yVariables = yVariables(:, cellfun(@isempty, regexpi(yVariables.Properties.VariableNames, 'mean$')));
yVariables = yVariables(:, cellfun(@isempty, regexpi(yVariables.Properties.VariableNames, 'meanlh$')));
yVariables = yVariables(:, cellfun(@isempty, regexpi(yVariables.Properties.VariableNames, 'meanrh$')));
ymat       = yVariables{:, :};
roiNames   = yVariables.Properties.VariableNames;
iid_concat = dataCThick.src_subject_id;
eid_concat = dataCThick.eventname;

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
             
% Show that rank deficiency exists even for MATLAB's automated modeling as
% it uses n-1 for all categorical variables
% tbl = cell2table([num2cell([dataWorking.smri_thick_cdk_fusiformrh, ageRecruitment, ageDelta]),           ...
%                  dataWorking.sex, imgVars.mri_info_deviceserialnumber, imgVars.mri_info_softwareversion, ...
%                  imgVars.household_income, imgVars.high_educ, num2cell(geneticPCs)], 'VariableNames',    ...
%                  [{'Thick', 'age', 'ageDelta', 'sex', 'device', 'software', 'income', 'education'}, strcat({'PC'}, num2str((1:20)', '%02d'))']);
% 
% mdl = fitlm(tbl, 'Thick ~ 1 + age + ageDelta + sex + device + software + income + education');

%% Create a new version of the design matrix with covariates
% First four columns should be: src_subject_id, eventname, rel_family_id, age
toIntersect = cell2table([designMatrix.src_subject_id, designMatrix.eventname, num2cell([designMatrix.rel_family_id, designMatrix.interview_age, covariates])], 'VariableNames', ...
                         [{'src_subject_id', 'eventname', 'rel_family_id', 'age'}, covarNames]);

% Save as T1w design matrix
writetable(toIntersect, fullfile(dirOut, 'DesignMatrix_T1w_ROI_CorticalThickness.csv'));

%% Intersect the design matrix and the data
[X, iid, eid, fid, agevec, ymatUse, contrasts, colnames_model, pihatmat, PregID, HomeID] = FEMA_intersect_design(fullfile(dirOut, 'DesignMatrix_T1w_ROI_CorticalThickness.csv'), ...
                                                                                                                 ymat, iid_concat, eid_concat);

%% Standardize ymatUse
ymatUse_std = (ymatUse - mean(ymatUse))./std(ymatUse);

%% Prepare random effects
SubjectEffect = iid;
FamilyEffect  = cellfun(@(x) strrep(x, ' ', ''), strcat({'F'}, num2str(fid)), 'UniformOutput', false);

%% Additional variables for FEMA_fit
nperms          = 100;
RandomEffects   = {'F', 'S', 'E'};
niter           = 1;
nbins           = 20;

%% Save variables
save(fullfile(dirOut, 'vars_analysis.mat'), 'SubjectEffect', 'FamilyEffect', 'X', 'iid', 'fid', 'agevec', 'ymatUse', 'ymatUse_std', 'colnames_model', '-v7.3');

%% Call FEMA_fit
initFEMA = tic;
[beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, Hessmat, logLikvec,                           ...
 beta_hat_perm, beta_se_perm, zmat_perm, sig2tvec_perm, sig2mat_perm, logLikvec_perm] =             ...
 FEMA_fit(covariates, SubjectEffect, eid, FamilyEffect, agevec, ymatUse_std, niter, contrasts,      ...
          nbins, [], 'RandomEffects', RandomEffects, 'nperms', nperms, 'PermType', 'wildbootstrap-nn');
elapsedFEMA = toc(initFEMA);

% Compare with MATLAB
mdl             = cell(size(ymatUse_std,2),1);
beta_MATLAB     = zeros(size(covariates,2),     size(ymatUse_std,2));
beta_se_MATLAB  = zeros(size(covariates,2),     size(ymatUse_std,2));
varComp_MATLAB  = zeros(length(RandomEffects),  size(ymatUse_std,2));
lowVar_MATLAB   = zeros(length(RandomEffects),  size(ymatUse_std,2));
UppVar_MATLAB   = zeros(length(RandomEffects),  size(ymatUse_std,2));
initMATLAB      = tic;
for mdls        = 1:size(ymatUse_std,2)
    mdl{mdls,1} = fitlmematrix(covariates, ymatUse_std(:,mdls), [{ones(length(FamilyEffect),1)}, {ones(length(SubjectEffect),1)}], {FamilyEffect, SubjectEffect});
end
elapsedMATLAB = toc(initMATLAB);

%% Extract parameters out of MATLAB
for mdls = 1:size(ymatUse_std,2)
    [a, b, c]               = covarianceParameters(mdl{mdls});
    beta_MATLAB(:, mdls)    = mdl{mdls}.Coefficients.Estimate;
    beta_se_MATLAB(:, mdls) = mdl{mdls}.Coefficients.SE;
    varComp_MATLAB(:, mdls) = [c{1}.Estimate^2, c{2}.Estimate^2, c{3}.Estimate^2];
    lowVar_MATLAB(:, mdls)  = [c{1}.Lower^2,    c{2}.Lower^2,    c{3}.Lower^2];
    UppVar_MATLAB(:, mdls)  = [c{1}.Upper^2,    c{2}.Upper^2,    c{3}.Upper^2];
end

%% Save everything
if not(exist(dirOut, 'dir'))
    mkdir(dirOut);
end
save(fullfile(dirOut, 'Results.mat'), '-v7.3');