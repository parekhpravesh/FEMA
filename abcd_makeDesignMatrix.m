%% Create design matrix - no filtering yet
%% Setup
% Specify release number
dataRelease = '4.0';

% Various paths
dirABCD         = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/';
dirCode         = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/cmig_tools-2.3.0';
dirOut          = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_re/home/pparekh/analyses/2023-02-17_FEMA-ABCD/CorticalThickness_vertexWise_FSE';
dirSupport      = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/support_files/ABCD_rel4.0_unfiltered/';
dirTabulated    = fullfile(dirABCD, dataRelease, 'tabulated', 'released'); 
fname_GRM       = fullfile(dirABCD, dataRelease, 'genomics', ['ABCD_rel', dataRelease, '_grm.mat']);

%% Read various files
% First, read the QC variables
dataQC = readtable(fullfile(dirABCD, dataRelease, 'tabulated', 'img', 'abcd_imgincl01.csv'));

% Next, get some additional info about MRI acquisition
dataInfo = readtable(fullfile(dirABCD, dataRelease, 'tabulated', 'img', 'abcd_mri01.csv'));

% Read various support files - only need allPCs and imgVars
imgVars       = readtable(fullfile(dirSupport, 'ABCD_rel4.0_covars_img_base_2yr.txt'));
allPCs        = readtable(fullfile(dirSupport, 'ABCD_rel4.0_pcs_base_2yr.txt'));

%% Make a merged field for aligning tables
dataQC.toAlign        = strcat(dataQC.src_subject_id,     dataQC.eventname);
imgVars.toAlign       = strcat(imgVars.src_subject_id,    imgVars.eventname);
allPCs.toAlign        = strcat(allPCs.src_subject_id,     allPCs.eventname);

%% Merge various support tables - imgVars and allPcs
locGenesisPCs = ~cellfun(@isempty, regexpi(allPCs.Properties.VariableNames, '^genesis'));
if sum(strcmpi(imgVars.toAlign, allPCs.toAlign)) ~= height(imgVars)
    error('Imaging and PC tables are not aligned');
end
imgVars = [imgVars, allPCs(:,locGenesisPCs)];

%% There are four subjects with two entries in QC file - get rid of them!
[a, b]   = histcounts(categorical(dataQC.toAlign));
toDelete = b(a == 2);
locs     = zeros(length(toDelete), 1);
for ids  = 1:length(toDelete)
    locs(ids,1) = find(strcmpi(dataQC.toAlign, toDelete{ids}), 1);
end
dataQC(locs, :) = [];

%% Align and merge dataQC table to imgVars
% Intersect tables
[~, b, c]  = intersect(imgVars.toAlign, dataQC.toAlign);

% Retain intersecting rows
imgVars = imgVars(b, :);
dataQC  = dataQC(c, :);

% Sanity check
if sum(strcmpi(dataQC.toAlign, imgVars.toAlign)) ~= height(dataQC)
    error('QC and imaging variable tables are not aligned');
end

% Merge tables
% At this stage, imgVars has the info on how subjects are ordered
locImgInclude = ~cellfun(@isempty, regexpi(dataQC.Properties.VariableNames, '^imginc'));
imgVars       = [imgVars, dataQC(:, locImgInclude)];

%% Subset data into baseline and follow up
% Only analyze subjects who have both the data points 
locWorking   = strcmpi(imgVars.eventname, 'baseline_year_1_arm_1');
dataBaseline = imgVars(locWorking,  :);
dataFollowup = imgVars(~locWorking, :);

% Only retain subjects who exist in both lists
[~, b, c]    = intersect(dataBaseline.src_subject_id, dataFollowup.src_subject_id);
dataBaseline = dataBaseline(b, :);
dataFollowup = dataFollowup(c, :);

% Ensure ordering matches
if sum(strcmpi(dataBaseline.src_subject_id, dataFollowup.src_subject_id)) ~= height(dataBaseline)
    error('Something went wrong during subsetting of data');
end

% Merge these into a working set
% At this stage dataWorking and yWorking have all the information we need
dataWorking = [dataBaseline; dataFollowup];

%% Get rid of missing information
% age at recruitment
% sex
% scanner device
% scanner software number
% household income
% highest parental education level
% first 20 genetic PCs (Genesis)
% Use dataWorking for all operations

% Location for PCs
locGenesisPCs = ~cellfun(@isempty, regexpi(dataWorking.Properties.VariableNames, '^genesis'));

% Unique event names - useful later
allUqEvents = unique(dataWorking.eventname);

% Find missing information
locMissing = isnan(dataWorking.interview_age)                                        | ...
             ismember(dataWorking.sex, {'NaN', 'NA'})                                | ...
             ismember(dataWorking.mri_info_deviceserialnumber, {'NaN', 'NA'})        | ...
             ismember(dataWorking.mri_info_softwareversion, {'NaN', 'NA'})           | ...
             ismember(dataWorking.household_income, {'NaN', 'NA'})                   | ...
             ismember(dataWorking.high_educ, {'NaN', 'NA'})                          | ...
             logical(sum(isnan(dataWorking{:, locGenesisPCs}),2));
missingIDs = dataWorking.src_subject_id(locMissing);

% Create toAlign variable from these missing IDs; loop over if more events
toDelete = [strcat(missingIDs, allUqEvents{1}); ...
            strcat(missingIDs, allUqEvents{2})];
        
% Delete these subjects from ymat and dataWorking
locDelete                 = ismember(dataWorking.toAlign, toDelete);
dataWorking(locDelete, :) = [];

%% Save
writetable(dataWorking, '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/GlobalDesignMatrix.csv', 'Delimiter', '\t');