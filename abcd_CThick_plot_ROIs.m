%% Read data
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/CorticalThickness_DK_FSE_nn';
data = load(fullfile(workDir, 'Summary_DK40.mat'));

% data    = readtable(fullfile(workDir, 'Summary_DK40.xlsx'));
% %% New roi names
% roiNames_Use = data.ROIName;

%% New names
roiNames = {'Bank STS', 'Caud Ant Cing', 'Caud Mid Fron', 'Cuneus', 'Entorhinal', 'Fusiform',   ...
            'Inf Parietal', 'Inf Temp', 'Isthmus Cing', 'Lat Occ', 'Lat Orb Fron', 'Lingual',  ...
            'Med Orb Fron', 'Mid Temp', 'Parahippocampal', 'Paracentral', 'Pars Opercularis',   ...
            'Pars orbitalis', 'Pars triangularis', 'Pericalcarine', 'Postcentral', 'Post Cing', ...
            'Precentral', 'Precuneus', 'Rost Ant Cing', 'Rost Mid Fron', 'Sup Fron',            ...
            'Sup Parietal', 'Sup Temp', 'Supramarginal', 'Fron Pole', 'Temp Pole', 'Trans Temp', 'Insula'};

%% Identify left and right hemispheres
locLeft     = find(~cellfun(@isempty, regexpi(data.roiNames_Use, 'lh$')));
locRight    = find(~cellfun(@isempty, regexpi(data.roiNames_Use, 'rh$')));

%% Plot all parameter estimates - age delta
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16], 'PaperUnits','centimeters', 'PaperSize', [16 16]);
allH = tight_subplot(1, 2, [0.0 0.02], [0.12 0.03], [0.195 0.01]);

% Color scheme
col_fitlme  = [217,95,2]./255;
col_FEMA    = [27,158,119]./255;

% Locations to plot at
xx1 = 0.7:1:34.7;
xx2 = 1.3:1:35.3;

% Hold everything
hold(allH(:), 'on');

% Settings
markerType  = 'o';
markerSize  = 3;
LineStyle   = 'none';
LineWidth   = 0.8;
capSize     = 4;

alpha = 0.05/68;

% MATLAB values
pValues_MATLAB_ageDelta  = data.pValues_MATLAB(data.locAgeDelta,:);
beta_MATLAB_ageDelta     = data.beta_MATLAB(data.locAgeDelta,:);
beta_se_MATLAB_ageDelta  = data.beta_se_MATLAB(data.locAgeDelta,:);

% FEMA values
pValues_FEMA_ageDelta    = data.pValues_FEMA(data.locAgeDelta,:);
beta_FEMA_ageDelta       = data.beta_hat(data.locAgeDelta,:);
beta_se_FEMA_ageDelta    = data.beta_se(data.locAgeDelta,:);

% Plot

% Slowly loop over all LH ROIs
for rois = 1:length(locLeft)    
    % Determine if statistically significant or not - FEMA
    if pValues_FEMA_ageDelta(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_FEMA_ageDelta(rois), xx2(rois), beta_se_FEMA_ageDelta(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_FEMA_ageDelta(rois), xx2(rois), beta_se_FEMA_ageDelta(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB_ageDelta(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_MATLAB_ageDelta(locLeft(rois)), xx1(rois), beta_se_MATLAB_ageDelta(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_MATLAB_ageDelta(locLeft(rois)), xx1(rois), beta_se_MATLAB_ageDelta(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
    end
end

% Slowly loop over all RH rois
for rois = 1:length(locRight)

    % Determine if statistically significant or not - FEMA
    if pValues_FEMA_ageDelta(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_FEMA_ageDelta(locRight(rois)), xx2(rois), beta_se_FEMA_ageDelta(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_FEMA_ageDelta(locRight(rois)), xx2(rois), beta_se_FEMA_ageDelta(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB_ageDelta(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_MATLAB_ageDelta(locRight(rois)), xx1(rois), beta_se_MATLAB_ageDelta(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_MATLAB_ageDelta(locRight(rois)), xx1(rois), beta_se_MATLAB_ageDelta(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
    end
end

for ax = 1:2
    % Apply axes settings
    allH(ax).YLim            = [0 35];
    allH(ax).Box             = 'off';
    allH(ax).XLim            = [-0.03 0.005];
    allH(ax).YTick           = 1:34;
    allH(ax).YGrid           = 'on';
    allH(ax).XLabel.String   = 'Estimated \beta age_{delta}';
    allH(ax).Title.FontWeight = 'normal';
    allH(ax).XAxis.FontSize  = 11;
    allH(ax).XLabel.FontSize = 12;
    allH(ax).YLabel.FontSize = 12;
    allH(ax).Title.FontSize  = 11;
    allH(ax).XTick           = -0.03:0.01:0.01;
    allH(ax).XTickLabel      = num2str((-0.03:0.01:0.01)', '%.2f');
end

allH(1).YTick           = 1:34;
allH(1).YTickLabel      = roiNames;
allH(1).YAxis.FontName  = 'CMUBright';
allH(1).XAxis.FontName  = 'CMUBright';
allH(2).YAxis.FontName  = 'CMUBright';
allH(2).XAxis.FontName  = 'CMUBright';

% Make a dashed zero line at zero
plot(allH(1), zeros(36), 0:35, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
plot(allH(2), zeros(36), 0:35, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);

allH(1).Title.String    = 'Left hemisphere';
allH(2).Title.String    = 'Right hemisphere';

allH(1).Title.FontName  = 'CMUBright';
allH(2).Title.FontName  = 'CMUBright';

% Overall legend
ll = legend(allH(2), {'FEMA', 'fitlmematrix'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off', 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2) - 0.12;
ll.Position(1) = allH(2).Position(1) - 0.25;

% Print figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/DK40_CThickness_ageDelta.png', '-dpng', '-r900');
close(fig);

% Clear up to prevent accidents!
clear *_ageDelta

%% Age of recruitment
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16], 'PaperUnits','centimeters', 'PaperSize', [16 16]);
allH = tight_subplot(1, 2, [0.0 0.02], [0.12 0.03], [0.195 0.01]);

% Color scheme
col_fitlme  = [217,95,2]./255;
col_FEMA    = [27,158,119]./255;

% Locations to plot at
xx1 = 0.7:1:34.7;
xx2 = 1.3:1:35.3;

% Hold everything
hold(allH(:), 'on');

% Settings
markerType  = 'o';
markerSize  = 3;
LineStyle   = 'none';
LineWidth   = 0.8;
capSize     = 4;

alpha = 0.05/68;

% MATLAB values
pValues_MATLAB_age  = data.pValues_MATLAB(data.locAge,:);
beta_MATLAB_age     = data.beta_MATLAB(data.locAge,:);
beta_se_MATLAB_age  = data.beta_se_MATLAB(data.locAge,:);

% FEMA values
pValues_FEMA_age    = data.pValues_FEMA(data.locAge,:);
beta_FEMA_age       = data.beta_hat(data.locAge,:);
beta_se_FEMA_age    = data.beta_se(data.locAge,:);

% Plot
% Slowly loop over all LH ROIs
for rois = 1:length(locLeft)    
    % Determine if statistically significant or not - FEMA
    if pValues_FEMA_age(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_FEMA_age(rois), xx2(rois), beta_se_FEMA_age(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_FEMA_age(rois), xx2(rois), beta_se_FEMA_age(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB_age(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_MATLAB_age(locLeft(rois)), xx1(rois), beta_se_MATLAB_age(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_MATLAB_age(locLeft(rois)), xx1(rois), beta_se_MATLAB_age(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
    end
end

% Slowly loop over all RH rois
for rois = 1:length(locRight)

    % Determine if statistically significant or not - FEMA
    if pValues_FEMA_age(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_FEMA_age(locRight(rois)), xx2(rois), beta_se_FEMA_age(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_FEMA_age(locRight(rois)), xx2(rois), beta_se_FEMA_age(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB_age(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_MATLAB_age(locRight(rois)), xx1(rois), beta_se_MATLAB_age(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_MATLAB_age(locRight(rois)), xx1(rois), beta_se_MATLAB_age(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
    end
end

for ax = 1:2
    % Apply axes settings
    allH(ax).YLim            = [0 35];
    allH(ax).Box             = 'off';
    allH(ax).XLim            = [-0.03 0.005];
    allH(ax).YTick           = 1:34;
    allH(ax).YGrid           = 'on';
    allH(ax).XLabel.String   = 'Estimated \beta age_{recruitment}';
    allH(ax).Title.FontWeight = 'normal';
    allH(ax).XAxis.FontSize  = 11;
    allH(ax).XLabel.FontSize = 12;
    allH(ax).YLabel.FontSize = 12;
    allH(ax).Title.FontSize  = 11;
    allH(ax).XTick           = -0.03:0.01:0.01;
    allH(ax).XTickLabel      = num2str((-0.03:0.01:0.01)', '%.2f');
end

allH(1).YTick           = 1:34;
allH(1).YTickLabel      = roiNames;
allH(1).YAxis.FontName  = 'CMUBright';
allH(1).XAxis.FontName  = 'CMUBright';
allH(2).YAxis.FontName  = 'CMUBright';
allH(2).XAxis.FontName  = 'CMUBright';

% Make a dashed zero line at zero
plot(allH(1), zeros(36), 0:35, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
plot(allH(2), zeros(36), 0:35, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);

allH(1).Title.String    = 'Left hemisphere';
allH(2).Title.String    = 'Right hemisphere';
allH(1).Title.FontName  = 'CMUBright';
allH(2).Title.FontName  = 'CMUBright';

% Overall legend
ll = legend(allH(2), {'FEMA', 'fitlmematrix'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off', 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2) - 0.12;
ll.Position(1) = allH(2).Position(1) - 0.25;

% Print figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/DK40_CThickness_ageRecruitment.png', '-dpng', '-r900');
close(fig);

% Clear up!
clear *_age

%% Random effects
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 22], 'PaperUnits','centimeters', 'PaperSize', [16 22]);
allH = tight_subplot(1, 3, [0.0 0.02], [0.09 0.03], [0.195 0.01]);

% Color scheme
col_fitlme  = [217,95,2]./255;
col_FEMA    = [27,158,119]./255;

% Locations to plot at
xx1 = 0.7:1:67.7;
xx2 = 1.3:1:68.3;

% Hold everything
hold(allH(:), 'on');

% Settings
markerType  = 'o';
markerSize  = 3;
LineStyle   = 'none';
LineWidth   = 0.8;
capSize     = 4;

% MATLAB values
% RFX_MATLAB      = [data.varComp_MATLAB(1,:); data.varComp_MATLAB(2,:); data.varComp_MATLAB(3,:)];
% RFX_low_MATLAB  = [data.lowFamily_MATLAB, data.lowSubject_MATLAB, data.lowError_MATLAB];
% RFX_upp_MATLAB  = [data.uppFamily_MATLAB, data.uppSubject_MATLAB, data.uppError_MATLAB];
% 
% % FEMA values
% RFX_FEMA        = [data.varFamily_FEMA, data.varSubject_FEMA, data.varError_FEMA];
% RFX_low_FEMA    = [data.lowFamily_FEMA, data.lowSubject_FEMA, data.lowError_FEMA];
% RFX_upp_FEMA    = [data.uppFamily_FEMA, data.uppSubject_FEMA, data.uppError_FEMA];

% Plot
for ax = 1:3
    errorbar(allH(ax), data.sig2mat(ax,:),        xx2, data.sig2Low(ax,:),       data.sig2Upp(ax,:),       'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    errorbar(allH(ax), data.varComp_MATLAB(ax,:), xx1, data.lowVar_MATLAB(ax,:), data.UppVar_MATLAB(ax,:), 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
end

allH(1).YTickLabel      = [strcat(roiNames, {' lh'}), strcat(roiNames, {' rh'})];
allH(1).YAxis.FontSize  = 8.5;
allH(1).YAxis.FontName  = 'CMUBright';
allH(1).Title.String    = 'Family effect (\itF\rm)';
allH(2).Title.String    = 'Subject effect (\itS\rm)';
allH(3).Title.String    = 'Unmodeled variance (\itE\rm)';
allH(1).Title.FontName  = 'CMUBright';
allH(2).Title.FontName  = 'CMUBright';
allH(3).Title.FontName  = 'CMUBright';

for ax = 1:3
    % Apply axes settings
    allH(ax).YLim            = [0 69];
    allH(ax).Box             = 'off';
    allH(ax).XLim            = [-0.2 1.8];
    allH(ax).YTick           = 1:68;
    allH(ax).YGrid           = 'on';
    allH(ax).XLabel.String   = 'Estimated \sigma^2';
    allH(ax).Title.FontWeight = 'normal';
    allH(ax).XAxis.FontSize  = 11;
    allH(ax).XLabel.FontSize = 12;
    allH(ax).YLabel.FontSize = 12;
    allH(ax).Title.FontSize  = 11;
    allH(ax).XTick           = 0.0:0.4:1.8;
    allH(ax).XTickLabels     = 0.0:0.4:1.8;
end

allH(1).YAxis.TickLength = [0 0];
allH(2).YAxis.TickLength = [0 0];
allH(3).YAxis.TickLength = [0 0];

% Overall legend
ll = legend(allH(2), {'FEMA', 'fitlmematrix'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off', 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2) - 0.09;

% Print figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/DK40_CThickness_RFX.png', '-dpng', '-r900');
close(fig);