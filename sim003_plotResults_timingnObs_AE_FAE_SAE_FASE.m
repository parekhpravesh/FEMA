%% Plot results of simulation 003: computational time as a function of number of observations
% This part plots results for AE, FAE, SAE, and FASE models
%% Load all results
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/fromTSD/';
load(fullfile(workDir, 'Results_experiment03.mat'));

%% Make a plot
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 11]);
allH = tight_subplot(2, 2, [0.2 0.085], [0.15 0.05], [0.08 0.022]);

% Color scheme and settings
col_FEMA    = [27,158,119]./255;
markerSize  = 10;

hold(allH(:), 'on');

% AE models
plot(allH(1), results.NumObservations, results.FEMA_AE,   'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
title(allH(1), '\itA\rm and \itE\rm model', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'CMUBright');
text(allH(1), -1800, 23, 'a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% FAE models
plot(allH(2), results.NumObservations, results.FEMA_FAE,   'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
title(allH(2), '\itF\rm, \itA\rm, and \itE\rm model', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'CMUBright');
text(allH(2), -1800, 23, 'b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% SAE models
plot(allH(3), results.NumObservations, results.FEMA_SAE,   'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
title(allH(3), '\itS\rm, \itA\rm, and \itE\rm model', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'CMUBright');
text(allH(3), -1800, 23, 'c)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% FASE models
plot(allH(4), results.NumObservations, results.FEMA_FASE,   'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
title(allH(4), '\itF\rm, \itA\rm, \itS\rm, and \itE\rm model', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'CMUBright');
text(allH(4), -1800, 23, 'd)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% Customize axes
for ax = 1:4
    allH(ax).XLim            = [settings.nObservations(1)-1000, settings.nObservations(end)+500];
    allH(ax).YLim            = [0 21];
    allH(ax).XTick           = settings.nObservations;
    allH(ax).XTickLabel      = num2str(settings.nObservations');
    allH(ax).XAxis.FontSize  = 8;
    allH(ax).YAxis.FontSize  = 8;
    allH(ax).XAxis.TickLabelRotation = 90;
    allH(ax).YTick           = 0:5:20;
    allH(ax).YTickLabel      = 0:5:20;
    allH(ax).XLabel.String   = 'Number of observations';
    allH(ax).YLabel.String   = 'Time elapsed (s)';
    allH(ax).XLabel.FontSize = 10;
    allH(ax).YLabel.FontSize = 10;
    allH(ax).XAxis.FontName  = 'CMUBright';
    allH(ax).YAxis.FontName  = 'CMUBright';
end

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/Results04-Comparison_ComputationalTime_AE-FAE-SAE-FASE.png', '-dpng', '-r900');
close(fig);