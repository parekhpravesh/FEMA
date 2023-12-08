%% Plot results of simulation 003: computational time as a function of number of y variables
%% Load all results
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/export_TSD/';
load(fullfile(workDir, 'Results_experiment03_nYvars.mat'));

%% Plot
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 9]);
allH = tight_subplot(1, 1, 0, [0.22 0.01], [0.1 0.025]);

% Color scheme and settings
col_fitlme      = [217,  95,   2]./255;
col_FEMA        = [27,  158, 119]./255;
col_fitlmePar   = [117, 112, 179]./255;
markerSize      = 10;

semilogy(allH(1), results.NumYVars, (results.FEMA_FSE),       'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
hold(allH(1), 'on');
semilogy(allH(1), results.NumYVars, (results.fitlme_FSE),     'LineStyle', '-', 'Color', col_fitlme, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
semilogy(allH(1), results.NumYVars, (results.fitlme_parFSE),  'LineStyle', '-', 'Color', col_fitlmePar, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);
box(allH(1), 'off');

% Customize axis
minValue                 = 0;
maxValue                 = (round(max(results.fitlme_FSE)+1000, 0));
ax                       = 1;
allH(1).YLim             = [1 5500];
allH(ax).XAxis.FontSize  = 10;
allH(ax).YAxis.FontSize  = 10;
allH(ax).XAxis.TickLabelRotation = 0;
allH(ax).YTick           = [1, 5, 10, 25:25:100, 250:250:1000, 2500:2500:5000];
allH(ax).XLabel.String   = 'Number of \ity\rm variables';
allH(ax).YLabel.String   = 'Time elapsed (s)';
allH(ax).XLabel.FontSize = 12;
allH(ax).YLabel.FontSize = 12;
allH(ax).XAxis.FontName  = 'CMUBright';
allH(ax).YAxis.FontName  = 'CMUBright';

ll = legend(allH(ax), {'FEMA', 'fitlmematrix', 'fitlmematrix (parallel)'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off', 'FontSize', 12, 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2) - 0.2;

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results03-Comparison_ComputationalTime_nYVars.png', '-dpng', '-r900');
close(fig);