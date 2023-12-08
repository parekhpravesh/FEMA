%% Reformat ROI information
load('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/Corrmat_FSE/tmp_roinames.mat');
gordonNames     = readtable('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/Corrmat_FSE/nodeNames.txt', 'ReadVariableNames', false);
gordonComm      = readtable('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/Corrmat_FSE/communityAffiliation.txt', 'ReadVariableNames', false);
gordonCommNames = readtable('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/Corrmat_FSE/communityNames.txt', 'ReadVariableNames', false);

% Parse Gordon parcellations
gordonNames.Var1      = strrep(strrep(gordonNames.Var1, 'L', 'Left '), 'R', 'Right ');
gordonNames.names     = strcat(gordonNames.Var1, gordonNames.Var2); 
gordonNames.Comm      = gordonComm.Var1;
gordonNames.CommNames = repmat({''}, height(gordonNames), 1);

for tmp = 1:height(gordonCommNames)
    loc = gordonNames.Comm == tmp;
    gordonNames.CommNames(loc) = gordonCommNames.Var1(tmp);
end

% Assign names to Gordon parcellations
roinames(217:549) = strcat(gordonNames.names, '-', gordonNames.CommNames);

% Sort Gordon names by community affiliation
[gordonNames, rc] = sortrows(gordonNames, 'Comm');

loc_LR_DK   = 1:68;
loc_LR_DT   = 69:216;
loc_LR_GP   = 217:549;
loc_LR_GP   = loc_LR_GP(rc);
loc_LR_aseg = 550:582;
loc_all     = [loc_LR_DK, loc_LR_DT, loc_LR_GP, loc_LR_aseg];

newnames = roinames(loc_all)';

%% Find locations
numParcels  = 582;
[r, c]      = find(tril(randn(numParcels, numParcels), -1));
idx         = sub2ind([numParcels numParcels], r, c);

%% Get Z scores
load('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/Corrmat_FSE/Results_half.mat', 'zmat');

%% Create names for atlas
atlasnames = [repmat({'DK40'}, length(loc_LR_DK), 1); repmat({'Destrieux'}, length(loc_LR_DT), 1); repmat({'Gordon'}, length(loc_LR_GP), 1); repmat({'aseg'}, length(loc_LR_aseg), 1)];

%% Put results together
results_age      = cell2table([atlasnames(r), atlasnames(c), newnames(r), newnames(c), num2cell(zmat(2,:)')], 'VariableNames', {'Atlas1', 'Atlas2', 'ROI1', 'ROI2', 'ZScore'});
results_ageDelta = cell2table([atlasnames(r), atlasnames(c), newnames(r), newnames(c), num2cell(zmat(3,:)')], 'VariableNames', {'Atlas1', 'Atlas2', 'ROI1', 'ROI2', 'ZScore'});

%% Write out
writetable(results_age,      '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/Corrmat_FSE/Lookup_Age.csv');
writetable(results_ageDelta, '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/2023-11-20_redone/Corrmat_FSE/Lookup_AgeDelta.csv');