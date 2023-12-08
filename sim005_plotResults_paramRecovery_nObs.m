%% Create MSE plots between FEMA and fitlmematrix
clear;
load('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/export_TSD/Results_experiment05.mat');

%% First plot - fixed effects
% Initialize figure
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
allH = tight_subplot(6, 1, [0.09 0.01], [0.15 0.04], [0.16 0.02]);

% Variable names
xText = {'X1', 'X2', 'X3', 'X4', 'X5', 'All fixed effects'};

% Average across repeats
avgTotalySqDiff_GTruth_FEMA_FFX = squeeze(mean(TotalySqDiff_GTruth_FEMA_FFX,2));
avgTotalySqDiff_GTruth_LME_FFX  = squeeze(mean(TotalySqDiff_GTruth_LME_FFX, 2));

% Actual plotting
for xvars = 1:5
    semilogx(allH(xvars), settings.nObservations, avgTotalySqDiff_GTruth_FEMA_FFX(xvars,:), 'LineStyle', '-', 'LineWidth', 2, 'Color', [27,158,119]./255,   'Marker', '.', 'MarkerSize', 10);
    hold(allH(xvars), 'on');
    semilogx(allH(xvars), settings.nObservations, avgTotalySqDiff_GTruth_LME_FFX(xvars,:),  'LineStyle', '-', 'LineWidth', 2, 'Color', [217,95,2,150]./255, 'Marker', '.', 'MarkerSize', 10);
end

% Add an overall plot - this is sum across variables
semilogx(allH(6), settings.nObservations, TotalmeanSqDiff_GTruth_FEMA_FFX, 'LineStyle', '-', 'LineWidth', 2, 'Color', [27,158,119]./255,   'Marker', '.', 'MarkerSize', 10);
hold(allH(6), 'on');
semilogx(allH(6), settings.nObservations, TotalmeanSqDiff_GTruth_LME_FFX, 'LineStyle', '-', 'LineWidth', 2, 'Color', [217,95,2,150]./255,   'Marker', '.', 'MarkerSize', 10);

% Customization
for xvars = 1:6
    % Turn off box
    box(allH(xvars), 'off');

    % Adjust X axis
    allH(xvars).XAxis.TickValues        = settings.nObservations;
    allH(xvars).XAxis.MinorTick         = 'off';
    allH(xvars).XAxis.TickLabelRotation = 90;
    allH(xvars).XAxis.FontSize          = 8;
    allH(xvars).XAxis.FontName          = 'CMUBright';

    % Some y axis settings
    allH(xvars).YAxis.FontSize          = 8;
    allH(xvars).YAxis.FontName          = 'CMUBright';
    allH(xvars).YAxis.TickLabelFormat   = '%.2f';
    
    % Add title
    title(allH(xvars), xText{xvars}, 'FontSize', 10, 'FontName', 'CMUBright', 'FontWeight', 'bold');
end

% Common Y axis
allH(3).YAxis.Label.String      = '$$\bf\sum_{y=1}^{500} (\beta_{true} - \beta_{estimated})^2$$';
allH(3).YAxis.Label.Interpreter = 'latex';
allH(3).YAxis.Label.FontSize    = 14;
allH(3).YAxis.Label.Position(2) = 0;

% Common X axis
allH(end).XAxis.Label.String      = 'Number of observations';
allH(end).XAxis.Label.FontSize    = 12;

% Common legend
ll = legend(allH(end), {'FEMA', 'fitlmematrix'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off', 'FontSize', 12, 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2) - 0.05;

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results05-FFX-SampleSizeParamRecovery.png', '-dpng', '-r900');
close(fig);

%% Second plot - random effects
% Initialize figure
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
allH = tight_subplot(4, 1, [0.09 0.01], [0.15 0.04], [0.16 0.02]);

% Variable names
xText = {'Family effect', 'Subject effect', 'Unmodeled variance', 'All random effects'};

% Average across repeats
avgTotalySqDiff_GTruth_FEMA_RFX = squeeze(mean(TotalySqDiff_GTruth_FEMA_RFX,2));
avgTotalySqDiff_GTruth_LME_RFX  = squeeze(mean(TotalySqDiff_GTruth_LME_RFX, 2));

% Actual plotting
for xvars = 1:3
    semilogx(allH(xvars), settings.nObservations, avgTotalySqDiff_GTruth_FEMA_RFX(xvars,:), 'LineStyle', '-', 'LineWidth', 2, 'Color', [27,158,119]./255,   'Marker', '.', 'MarkerSize', 10);
    hold(allH(xvars), 'on');
    semilogx(allH(xvars), settings.nObservations, avgTotalySqDiff_GTruth_LME_RFX(xvars,:),  'LineStyle', '-', 'LineWidth', 2, 'Color', [217,95,2,150]./255, 'Marker', '.', 'MarkerSize', 10);
end

% Add an overall plot - this is sum across variables
semilogx(allH(4), settings.nObservations, TotalmeanSqDiff_GTruth_FEMA_RFX, 'LineStyle', '-', 'LineWidth', 2, 'Color', [27,158,119]./255,   'Marker', '.', 'MarkerSize', 10);
hold(allH(4), 'on');
semilogx(allH(4), settings.nObservations, TotalmeanSqDiff_GTruth_LME_RFX, 'LineStyle', '-', 'LineWidth', 2, 'Color', [217,95,2,150]./255,   'Marker', '.', 'MarkerSize', 10);

% Customization
for xvars = 1:4
    % Turn off box
    box(allH(xvars), 'off');

    % Adjust X axis
    allH(xvars).XAxis.TickValues        = settings.nObservations;
    allH(xvars).XAxis.MinorTick         = 'off';
    allH(xvars).XAxis.TickLabelRotation = 90;
    allH(xvars).XAxis.FontSize          = 8;
    allH(xvars).XAxis.FontName          = 'CMUBright';

    % Some y axis settings
    allH(xvars).YAxis.FontSize          = 8;
    allH(xvars).YAxis.FontName          = 'CMUBright';
    allH(xvars).YAxis.TickLabelFormat   = '%.2f';
    
    % Add title
    title(allH(xvars), xText{xvars}, 'FontSize', 10, 'FontName', 'CMUBright', 'FontWeight', 'bold');
end

% Common Y axis
allH(2).YAxis.Label.String      = '$$\bf\sum_{y=1}^{500} (\sigma_{true}^2 - \sigma_{estimated}^2)^2$$';
allH(2).YAxis.Label.Interpreter = 'latex';
allH(2).YAxis.Label.FontSize    = 14;
allH(2).YAxis.Label.Position(2) = 0;

% Common X axis
allH(end).XAxis.Label.String      = 'Number of observations';
allH(end).XAxis.Label.FontSize    = 12;

% Common legend
ll = legend(allH(4), {'FEMA', 'fitlmematrix'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off', 'FontSize', 12, 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2) - 0.12;

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results05-RFX-SampleSizeParamRecovery.png', '-dpng', '-r900');
close(fig);

%% Deprecated code for differences
% %% Initialize figure
% fig  = figure('Units', 'centimeters', 'Position', [10 10 16 8]);
% allH = tight_subplot(1, 1, 0, [0.23 0.02], [0.16 0.02]);
% % allH = tight_subplot(1, 1, [0.17, 0.15], [0.12 0.03], [0.16 0.02]);
% 
% % First plot - fixed effects
% for xvars = 1:5
%     semilogx(allH(1), settings.nObservations, totalSqDiff_FEMA_LME(xvars, :), 'LineStyle', ':', 'LineWidth', 1, 'Color', [0 0 0], 'Marker', '.');
%     hold(allH(1), 'on');
% end
% % Average
% semilogx(allH(1), settings.nObservations, AvgTotalSqDiff_FEMA_LME, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [217 95 2]./255, 'Marker', '.', 'MarkerSize', 14);
% 
% % Customize axis
% box(allH(1), 'off');
% 
% % Adjust X axis
% allH(1).XAxis.TickValues        = settings.nObservations;
% allH(1).XAxis.MinorTick         = 'off';
% allH(1).XAxis.TickLabelRotation = 90;
% allH(1).XAxis.FontSize          = 10;
% allH(1).XAxis.Label.String      = 'Number of observations';
% allH(1).XAxis.Label.FontSize    = 12;
% allH(1).XAxis.FontName          = 'CMUBright';
% 
% % Adjust Y axis
% allH(1).YAxis.FontSize          = 10;
% allH(1).YAxis.Label.String      = '$$\bf\sum_{y=1}^{500} (\beta_{FEMA} - \beta_{fitlmematrix})^2$$';
% allH(1).YAxis.Label.Interpreter = 'latex';
% allH(1).YAxis.Label.FontSize    = 12;
% allH(1).YAxis.TickLabelFormat   = '%.2f';
% allH(1).YAxis.FontName          = 'CMUBright';
% 
% % Adjust title
% % allH(1).Title.String     = 'Fixed Effects';
% % allH(1).Title.FontName   = 'CMUBright';
% 
% print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results05-FFX-SampleSizeParamRecovery.png', '-dpng', '-r900');
% close(fig);
% 
% %% Second plot - random effects
% fig  = figure('Units', 'centimeters', 'Position', [10 10 16 8]);
% allH = tight_subplot(1, 1, 0, [0.23 0.02], [0.16 0.02]);
% 
% for xvars = 1:3
%     semilogx(allH(1), settings.nObservations, totalSqDiff_RFX_FEMA_LME(xvars, :), 'LineStyle', ':', 'LineWidth', 1, 'Color', [0 0 0], 'Marker', '.');
%     hold(allH(1), 'on');
% end
% 
% % Average
% semilogx(allH(1), settings.nObservations, AvgTotalSqDiff_RFX_FEMA_LME, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [217 95 2]./255, 'Marker', '.', 'MarkerSize', 14);
% 
% % Customize axis
% box(allH(1), 'off');
% 
% % Adjust X axis
% allH(1).XAxis.TickValues        = settings.nObservations;
% allH(1).XAxis.MinorTick         = 'off';
% allH(1).XAxis.TickLabelRotation = 90;
% allH(1).XAxis.FontSize          = 10;
% allH(1).XAxis.Label.String      = 'Number of observations';
% allH(1).XAxis.Label.FontSize    = 12;
% allH(1).XAxis.FontName          = 'CMUBright';
% 
% % Adjust Y axis
% allH(1).YAxis.FontSize          = 10;
% allH(1).YAxis.Label.String      = '$$\bf\sum_{y=1}^{500} (\sigma_{FEMA}^2 - \sigma_{fitlmematrix}^2)^2$$';
% allH(1).YAxis.Label.Interpreter = 'latex';
% allH(1).YAxis.Label.FontSize    = 12;
% allH(1).YAxis.TickLabelFormat   = '%.2f';
% allH(1).XAxis.FontName          = 'CMUBright';
% 
% % Adjust title
% % allH(2).Title.String    = 'Random Effects';
% % allH(2).Title.FontName  = 'CMUBright';
% 
% % Put a) and b)
% % text(allH(1), 25, allH(1).Title.Position(2),  'a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');
% % text(allH(2), 25, allH(2).Title.Position(2),  'b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');
% 
% %% Save figure
% print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results05-RFX-SampleSizeParamRecovery.png', '-dpng', '-r900');
% close(fig);