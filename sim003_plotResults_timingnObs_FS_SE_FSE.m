%% Plot results of simulation 003: computational time as a function of number of observations
% This part plots results for FS, SE, and FSE models, comparing with
% fitlmematrix timing
%% Load all results
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/export_TSD/';
load(fullfile(workDir, 'Results_experiment03_nObs.mat'));

%% Make a plot
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 9], 'PaperUnits','centimeters', 'PaperSize',[16 16]);
allH = tight_subplot(1, 3, [0.105 0.085], [0.24 0.06], [0.08 0.022]);

% Color scheme and settings
col_fitlme      = [217, 95,  2]./255;
col_fitlmePar   = [117, 112, 179]./255;
col_FEMA        = [27,  158, 119]./255;
markerSize      = 10;

hold(allH(:), 'on');

% First figure - FE, SE, FSE
% FE models
plot(allH(1), results.NumObservations, results.FEMA_FE,       'LineStyle', '-', 'Color', col_FEMA,      'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
plot(allH(1), results.NumObservations, results.fitlme_FE,     'LineStyle', '-', 'Color', col_fitlme,    'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
plot(allH(1), results.NumObservations, results.fitlme_par_FE, 'LineStyle', '-', 'Color', col_fitlmePar, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
title(allH(1), '\itF\rm and \itE\rm model', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'CMUBright');
text(allH(1), -4000, 75, 'a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% SE models
plot(allH(2), results.NumObservations, results.FEMA_SE,         'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
plot(allH(2), results.NumObservations, results.fitlme_SE,       'LineStyle', '-', 'Color', col_fitlme, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
plot(allH(2), results.NumObservations, results.fitlme_par_SE,   'LineStyle', '-', 'Color', col_fitlmePar, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
title(allH(2), '\itS\rm and \itE\rm model', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'CMUBright');
text(allH(2), -4000, 75, 'b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% FSE models
plot(allH(3), results.NumObservations, results.FEMA_FSE,        'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
plot(allH(3), results.NumObservations, results.fitlme_FSE,      'LineStyle', '-', 'Color', col_fitlme, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
plot(allH(3), results.NumObservations, results.fitlme_par_FSE,  'LineStyle', '-', 'Color', col_fitlmePar, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
title(allH(3), '\itF\rm, \itS\rm and \itE\rm model', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'CMUBright');
text(allH(3), -4000, 75, 'c)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% Customize axes
for ax = 1:3
    allH(ax).XLim            = [settings.nObservations(1)-1000, settings.nObservations(end)+500];
    allH(ax).YLim            = [0 72];
    allH(ax).XTick           = settings.nObservations;
    allH(ax).XTickLabel      = num2str(settings.nObservations');
    allH(ax).XAxis.FontSize  = 8;
    allH(ax).YAxis.FontSize  = 8;
    allH(ax).XAxis.TickLabelRotation = 90;
    allH(ax).YTick           = 0:10:70;
    allH(ax).YTickLabel      = 0:10:70;
    allH(ax).XLabel.String   = 'Number of observations';
    allH(ax).YLabel.String   = 'Time elapsed (s)';
    allH(ax).XLabel.FontSize = 10;
    allH(ax).YLabel.FontSize = 10;
    allH(ax).XAxis.FontName  = 'CMUBright';
    allH(ax).YAxis.FontName  = 'CMUBright';
end

% Add a legend
ll             = legend(allH(2), {'FEMA', 'fitlmematrix', 'fitlmematrix (parallel)'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off', 'FontSize', 12, 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2)-0.25;
ll.Position(1) = 0.15;

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results03-Comparison_ComputationalTime_FS-SE-FSE.png', '-dpng', '-r900');
close(fig);