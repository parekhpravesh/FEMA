%% Analyze corrmat
%% Setup
% Add FEMA code to path
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta/cmig_tools_utils/matlab');
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta/FEMA');

% Specify release number
dataRelease = '4.0';

% Various paths
dirABCD         = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/';
dirCode         = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta';
dirOut          = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE';
dirSupport      = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered/';
dirDesign       = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD';
dirTabulated    = fullfile(dirABCD, dataRelease, 'tabulated', 'released'); 
dirImaging      = fullfile(dirABCD, dataRelease, 'imaging_concat', 'corrmat', 'restingstate');
fnameGRM        = fullfile(dirABCD, dataRelease, 'genomics', ['ABCD_rel', dataRelease, '_grm.mat']);
fnameDesign     = 'GlobalDesignMatrix.csv';
fnameImaging    = 'rsfmri_fd0.20';

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

%% Get corrmat data
% [ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat] = FEMA_process_data(fnameImaging, dirTabulated, dirImaging, 'corrmat', 'pihat_file', fnameGRM);
[ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat] = parseCorrmat(fnameImaging, dirTabulated, dirImaging, fnameGRM);

% Get rid of subjects with one visit
[a, b]      = histcounts(categorical(iid_concat));
toDelete    = b(a == 1);
locs        = ismember(iid_concat, toDelete);
ymat        = ymat(~locs, :);
iid_concat  = iid_concat(~locs, :);
eid_concat  = eid_concat(~locs, :);

% Add a column for alignment
ymat_toAlign = strcat(iid_concat, eid_concat);

%% Within the design matrix, remove subjects who did not pass QC
toDelete = unique(designMatrix.src_subject_id(designMatrix.imgincl_rsfmri_include == 0));
designMatrix(ismember(designMatrix.src_subject_id, toDelete), :) = [];

% %% There are two subjects who are in designMatrix but not in iid_concat
% missingID = setdiff(designMatrix.src_subject_id, iid_concat);
% designMatrix(ismember(designMatrix.src_subject_id, missingID), :) = [];
% 
% Ensure that data lines up with design matrix
% [a, ~]       = ismember(ymat_toAlign, designMatrix.toAlign);
% ymat         = ymat(a, :);
% iid_concat   = iid_concat(a, :);
% eid_concat   = eid_concat(a, :);
% ymat_toAlign = ymat_toAlign(a, :);
% 
% [a, b]       = ismember(designMatrix.toAlign, ymat_toAlign);
% ymat         = ymat(b, :);
% iid_concat   = iid_concat(b, :);
% eid_concat   = eid_concat(b, :);
% ymat_toAlign = ymat_toAlign(b, :);

%% Intersect designMatrix and ymat_toAlign
[a, b, c]    = intersect(designMatrix.toAlign, ymat_toAlign);
designMatrix = designMatrix(b, :);
ymat_toAlign = ymat_toAlign(c, :);
ymat         = ymat(c, :);
iid_concat   = iid_concat(c, :);
eid_concat   = eid_concat(c, :);

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

% Save as design matrix
writetable(toIntersect, fullfile(dirOut, 'DesignMatrix_rsfMRI_corrmat_half.csv'));

%% Intersect the design matrix and the corrmat data
[X, iid, eid, fid, agevec, ymatUse, contrasts, colnames_model, pihatmat, PregID, HomeID] = FEMA_intersect_design(fullfile(dirOut, 'DesignMatrix_rsfMRI_corrmat_half.csv'), ...
                                                                                                              ymat, iid_concat, eid_concat, 'pihat', pihat);

%% Standardize ymatUse
ymatUse_std = (ymatUse - mean(ymatUse))./std(ymatUse);

%% Prepare random effects
SubjectEffect = iid;
FamilyEffect  = cellfun(@(x) strrep(x, ' ', ''), strcat({'F'}, num2str(fid)), 'UniformOutput', false);

%% Additional variables for FEMA_fit
nperms          = 0;
RandomEffects   = {'F', 'S', 'E'};
niter           = 1;
nbins           = 20;

%% Call FEMA_fit
initFEMA = tic;
[beta_hat, beta_se, zmat, logpmat, sig2tvec, sig2mat, Hessmat, logLikvec] =   ...
 FEMA_fit(X, SubjectEffect, eid, FamilyEffect, agevec, ymatUse_std, niter,    ...
          contrasts, nbins, [], 'RandomEffects', RandomEffects, 'nperms', nperms);
elapsedFEMA = toc(initFEMA);    

%% Save everything
if not(exist(dirOut, 'dir'))
    mkdir(dirOut);
end
save(fullfile(dirOut, 'Results_half.mat'), '-v7.3');

function [ymat, iid_concat, eid_concat, ivec_mask, mask, colnames_imaging, pihat] = parseCorrmat(fnameImaging, dirTabulated, dirImaging, fnameGRM)
% Code borrowed from FEMA_process_data
% Load data - takes forever
fname_corrmat   = sprintf('%s/%s.mat', dirImaging, fnameImaging);
tmp             = load(fname_corrmat);
nframes_min     = 375; % This should be optional  input param
ivec_tmp        = find(tmp.nsumvec>=nframes_min);

% Unset
colnames_imaging = [];

% Initialize a dummy matrix
nparcels = size(tmp.measmat,3);
[r, c]   = find(tril(randn(nparcels, nparcels), -1));
idx      = sub2ind([nparcels, nparcels], r, c);

% Create y variable
ymat      = tmp.measmat(ivec_tmp, idx);
dirlist   = tmp.dirlist(ivec_tmp);

% Initialize some variables
subjidvec   = cell(size(dirlist)); 
sitevec     = cell(size(dirlist)); 
datevec     = cell(size(dirlist)); 
visitidvec  = cell(size(dirlist));

% Populate these variables
for diri = 1:length(dirlist)
    tmp = regexp(dirlist{diri}, '^BOLDPROC_(?<site>\w+)_(?<SubjID>\w+)_(?<event>[^_]*)_(?<date>[^_.]+).(?<time>[^_]+)_', 'names');
    subjidvec{diri}  = tmp.SubjID;
    sitevec{diri}    = tmp.site;
    eventvec{diri}   = tmp.event;
    datevec{diri}    = tmp.date;
    visitidvec{diri} = sprintf('%s_%s_%s',tmp.site,tmp.SubjID,tmp.event);
end
iid_concat      = colvec(strcat('NDAR_',subjidvec)); % Not sure why imaging subject IDs are missing the NDAR_ part
visitid_concat  = visitidvec;
ivec_mask       = [];
mask            = [];

files           = dir([dirTabulated '/abcd_mri01*']);
fname_tabulated = [dirTabulated '/' files.name];

logging('Reading tabulated imaging data from %s',fname_tabulated);
imgtable         = readtable(fname_tabulated);
imgtable.idevent = strcat(imgtable.src_subject_id,'_',imgtable.eventname);

files            = dir([dirTabulated '/abcd_imgincl01*']);
fname_incflag    = [dirTabulated '/' files.name];
inctable         = readtable(fname_incflag);
inctable.idevent = strcat(inctable.src_subject_id,'_',inctable.eventname);

if ismember('dataset_id', imgtable.Properties.VariableNames)
    imgtable = removevars(imgtable,'dataset_id');
    inctable = removevars(inctable,'dataset_id');
end

[dummy IA IB]    = intersect(imgtable.idevent,inctable.idevent,'stable');
imgtable         = join(imgtable(IA,:),inctable(IB,:));

iid_imgtable     = imgtable.subjectkey;
eid_imgtable     = imgtable.eventname;
date_imgtable    = imgtable.mri_info_studydate;
visitid_imgtable = imgtable.mri_info_visitid;

% Merge vertex/voxelwise and tabulated imaging data
tmplist_concat      = visitid_concat;
tmplist_imgtable    = visitid_imgtable;
[dummy IA IB]       = intersect(tmplist_concat,tmplist_imgtable,'stable');
eid_concat          = eid_imgtable(IB);
iid_concat          = iid_concat(IA);
idevent             = strcat(iid_concat(:),'_',eid_concat(:));
imgtable            = imgtable(IB,:);
ymat                = ymat(IA,:);

% Get pihat
pihat = load(fnameGRM);
end