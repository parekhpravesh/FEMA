%% Create MSE plots between FEMA and fitlmematrix
clear;
load('/Applications/Projects/2023-01-31_FEMAExperiments/fromTSD/2023-08-02_redone/Summary_Simulation004_Compare_fitlme_perms_nn_repeats.mat');

%% Initialize figure
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 8]);
allH = tight_subplot(1, 1, 0, [0.23 0.02], [0.16 0.02]);
% allH = tight_subplot(1, 1, [0.17, 0.15], [0.12 0.03], [0.16 0.02]);

% First plot - fixed effects
for xvars = 1:5
    semilogx(allH(1), settings.nObservations, totalSqDiff_FEMA_LME(xvars, :), 'LineStyle', ':', 'LineWidth', 1, 'Color', [0 0 0], 'Marker', '.');
    hold(allH(1), 'on');
end
% Average
semilogx(allH(1), settings.nObservations, AvgTotalSqDiff_FEMA_LME, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [217 95 2]./255, 'Marker', '.', 'MarkerSize', 14);

% Customize axis
box(allH(1), 'off');

% Adjust X axis
allH(1).XAxis.TickValues        = settings.nObservations;
allH(1).XAxis.MinorTick         = 'off';
allH(1).XAxis.TickLabelRotation = 90;
allH(1).XAxis.FontSize          = 10;
allH(1).XAxis.Label.String      = 'Number of observations';
allH(1).XAxis.Label.FontSize    = 12;
allH(1).XAxis.FontName          = 'CMUBright';

% Adjust Y axis
allH(1).YAxis.FontSize          = 10;
allH(1).YAxis.Label.String      = '$$\bf\sum_{y=1}^{500} (\beta_{FEMA} - \beta_{fitlmematrix})^2$$';
allH(1).YAxis.Label.Interpreter = 'latex';
allH(1).YAxis.Label.FontSize    = 12;
allH(1).YAxis.TickLabelFormat   = '%.2f';
allH(1).YAxis.FontName          = 'CMUBright';

% Adjust title
% allH(1).Title.String     = 'Fixed Effects';
% allH(1).Title.FontName   = 'CMUBright';

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/Results06-FFX-SampleSizeParamRecovery.png', '-dpng', '-r900');
close(fig);

%% Second plot - random effects
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 8]);
allH = tight_subplot(1, 1, 0, [0.23 0.02], [0.16 0.02]);

for xvars = 1:3
    semilogx(allH(1), settings.nObservations, totalSqDiff_RFX_FEMA_LME(xvars, :), 'LineStyle', ':', 'LineWidth', 1, 'Color', [0 0 0], 'Marker', '.');
    hold(allH(1), 'on');
end
    
% Average
semilogx(allH(1), settings.nObservations, AvgTotalSqDiff_RFX_FEMA_LME, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [217 95 2]./255, 'Marker', '.', 'MarkerSize', 14);

% Customize axis
box(allH(1), 'off');

% Adjust X axis
allH(1).XAxis.TickValues        = settings.nObservations;
allH(1).XAxis.MinorTick         = 'off';
allH(1).XAxis.TickLabelRotation = 90;
allH(1).XAxis.FontSize          = 10;
allH(1).XAxis.Label.String      = 'Number of observations';
allH(1).XAxis.Label.FontSize    = 12;
allH(1).XAxis.FontName          = 'CMUBright';

% Adjust Y axis
allH(1).YAxis.FontSize          = 10;
allH(1).YAxis.Label.String      = '$$\bf\sum_{y=1}^{500} (\sigma_{FEMA}^2 - \sigma_{fitlmematrix}^2)^2$$';
allH(1).YAxis.Label.Interpreter = 'latex';
allH(1).YAxis.Label.FontSize    = 12;
allH(1).YAxis.TickLabelFormat   = '%.2f';
allH(1).XAxis.FontName          = 'CMUBright';

% Adjust title
% allH(2).Title.String    = 'Random Effects';
% allH(2).Title.FontName  = 'CMUBright';

% Put a) and b)
% text(allH(1), 25, allH(1).Title.Position(2),  'a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');
% text(allH(2), 25, allH(2).Title.Position(2),  'b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

%% Save figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/Results06-RFX-SampleSizeParamRecovery.png', '-dpng', '-r900');
close(fig);