%% Read data
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/fromCMIG/CorticalThickness_DK_FSE_nn';
data    = readtable(fullfile(workDir, 'Summary_DK40.xlsx'));

%% New roi names
roiNames_Use = data.ROIName;

%% New names
roiNames = {'Bank STS', 'Caud Ant Cing', 'Caud Mid Fron', 'Cuneus', 'Entorhinal', 'Fusiform',   ...
            'Inf Parietal', 'Inf Temp', 'Isthmus Cing', 'Lat Occ', 'Lat Orb Fron', 'Lingual',  ...
            'Med Orb Fron', 'Mid Temp', 'Parahippocampal', 'Paracentral', 'Pars Opercularis',   ...
            'Pars orbitalis', 'Pars triangularis', 'Pericalcarine', 'Postcentral', 'Post Cing', ...
            'Precentral', 'Precuneus', 'Rost Ant Cing', 'Rost Mid Fron', 'Sup Fron',            ...
            'Sup Parietal', 'Sup Temp', 'Supramarginal', 'Fron Pole', 'Temp Pole', 'Trans Temp', 'Insula'};

%% Identify left and right hemispheres
locLeft     = find(~cellfun(@isempty, regexpi(roiNames_Use, 'lh$')));
locRight    = find(~cellfun(@isempty, regexpi(roiNames_Use, 'rh$')));

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
pValues_MATLAB  = data.pValue_ageDelta_MATLAB;
beta_MATLAB     = data.Beta_ageDelta_MATLAB;
beta_se_MATLAB  = data.SE_ageDelta_MATLAB;

% FEMA values
pValues_FEMA    = data.pValue_ageDelta_FEMA;
beta_FEMA       = data.Beta_ageDelta_FEMA;
beta_se_FEMA    = data.SE_ageDelta_FEMA;

% Plot

% Slowly loop over all LH ROIs
for rois = 1:length(locLeft)    
    % Determine if statistically significant or not - FEMA
    if pValues_FEMA(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_FEMA(rois), xx2(rois), beta_se_FEMA(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_FEMA(rois), xx2(rois), beta_se_FEMA(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_MATLAB(locLeft(rois)), xx1(rois), beta_se_MATLAB(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_MATLAB(locLeft(rois)), xx1(rois), beta_se_MATLAB(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
    end
end

% Slowly loop over all RH rois
for rois = 1:length(locRight)

    % Determine if statistically significant or not - FEMA
    if pValues_FEMA(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_FEMA(locRight(rois)), xx2(rois), beta_se_FEMA(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_FEMA(locRight(rois)), xx2(rois), beta_se_FEMA(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_MATLAB(locRight(rois)), xx1(rois), beta_se_MATLAB(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_MATLAB(locRight(rois)), xx1(rois), beta_se_MATLAB(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
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
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/DK40_CThickness_ageDelta.png', '-dpng', '-r900');
close(fig);

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
pValues_MATLAB  = data.pValue_age_MATLAB;
beta_MATLAB     = data.Beta_age_MATLAB;
beta_se_MATLAB  = data.SE_age_MATLAB;

% FEMA values
pValues_FEMA    = data.pValue_age_FEMA;
beta_FEMA       = data.Beta_age_FEMA;
beta_se_FEMA    = data.SE_age_FEMA;

% Plot
% Slowly loop over all LH ROIs
for rois = 1:length(locLeft)    
    % Determine if statistically significant or not - FEMA
    if pValues_FEMA(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_FEMA(rois), xx2(rois), beta_se_FEMA(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_FEMA(rois), xx2(rois), beta_se_FEMA(rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB(locLeft(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(1), beta_MATLAB(locLeft(rois)), xx1(rois), beta_se_MATLAB(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(1), beta_MATLAB(locLeft(rois)), xx1(rois), beta_se_MATLAB(locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
    end
end

% Slowly loop over all RH rois
for rois = 1:length(locRight)

    % Determine if statistically significant or not - FEMA
    if pValues_FEMA(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_FEMA(locRight(rois)), xx2(rois), beta_se_FEMA(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_FEMA(locRight(rois)), xx2(rois), beta_se_FEMA(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
    end

    % Determine if statistically significant or not - MATLAB
    if pValues_MATLAB(locRight(rois)) < alpha
        
        % Plot with filled circle
        errorbar(allH(2), beta_MATLAB(locRight(rois)), xx1(rois), beta_se_MATLAB(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
    else
        % Plot without filling
        errorbar(allH(2), beta_MATLAB(locRight(rois)), xx1(rois), beta_se_MATLAB(locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
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
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/DK40_CThickness_ageRecruitment.png', '-dpng', '-r900');
close(fig);

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
RFX_MATLAB      = [data.varFamily_MATLAB, data.varSubject_MATLAB, data.varError_MATLAB];
RFX_low_MATLAB  = [data.lowFamily_MATLAB, data.lowSubject_MATLAB, data.lowError_MATLAB];
RFX_upp_MATLAB  = [data.uppFamily_MATLAB, data.uppSubject_MATLAB, data.uppError_MATLAB];

% FEMA values
RFX_FEMA        = [data.varFamily_FEMA, data.varSubject_FEMA, data.varError_FEMA];
RFX_low_FEMA    = [data.lowFamily_FEMA, data.lowSubject_FEMA, data.lowError_FEMA];
RFX_upp_FEMA    = [data.uppFamily_FEMA, data.uppSubject_FEMA, data.uppError_FEMA];

% Plot
for ax = 1:3
    errorbar(allH(ax), RFX_FEMA(:, ax),   xx2, RFX_low_FEMA(:, ax),   RFX_upp_FEMA(:, ax),   'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    errorbar(allH(ax), RFX_MATLAB(:, ax), xx1, RFX_low_MATLAB(:, ax), RFX_upp_MATLAB(:, ax), 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
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
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/DK40_CThickness_RFX.png', '-dpng', '-r900');
close(fig);