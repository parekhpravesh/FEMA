%% Generate figures for graphical abstract
%% As a function of number of y variables
% Load all results
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/export_TSD';
outDir  = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone';
load(fullfile(workDir, 'Results_experiment03_nYvars.mat'));

% Make a plot
fig  = figure('Units', 'centimeters', 'Position', [10 10 5.5 3.5], 'Color', [0 0 0], 'InvertHardcopy','off');
allH = tight_subplot(1, 1, [0 0], [0.34 0.04], [0.13 0.025]);

% Color scheme and settings
col_FEMA        = [27,  158, 119]./255;

hold(allH(:), 'on');

plot(allH(1), results.NumYVars, results.FEMA_FSE,   'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 5);

% Customize axis
ax                       = 1;
allH(ax).XLim            = [90 5010];
allH(ax).YLim            = [0 round(max(results.FEMA_FSE)+1, 0)];
allH(ax).XTick           = 500:500:5000;
allH(ax).XTickLabel      = 500:500:5000;
allH(ax).XAxis.FontSize  = 8;
allH(ax).YAxis.FontSize  = 8;
allH(ax).XAxis.TickLabelRotation = 90;
allH(ax).YTick           = 1:1:5;
allH(ax).YTickLabel      = 1:1:5;
% allH(ax).YTick           = 1:2:round(max(results.FEMA_FSE)+1, 0);
% allH(ax).YTickLabel      = 1:2:round(max(results.FEMA_FSE)+1, 0);
allH(ax).XLabel.String   = 'Number of \ity\rm variables';
allH(ax).YLabel.String   = 'Time elapsed (s)';
allH(ax).XLabel.FontSize = 9;
allH(ax).YLabel.FontSize = 9;
allH(ax).XColor = [1 1 1];
allH(ax).YColor = [1 1 1];
allH(ax).Color  = [0 0 0];

print(fullfile(outDir, 'GraphicalAbstract_numYvars.png'), '-dpng', '-r900');
close(fig);

%% Number of observations - FSE model
%% Load all results
load(fullfile(workDir, 'Results_experiment03_nObs.mat'));

% Make a plot
fig  = figure('Units', 'centimeters', 'Position', [10 10 5.5 3.5], 'Color', [0 0 0], 'InvertHardcopy','off');
allH = tight_subplot(1, 1, [0 0], [0.37 0.04], [0.13 0.025]);

% Color scheme and settings
col_FEMA        = [27,  158, 119]./255;
markerSize      = 5;

hold(allH(:), 'on');

% FSE models
plot(allH(1), results.NumObservations, results.FEMA_FSE, 'LineStyle', '-', 'Color', col_FEMA, 'LineWidth', 1, 'Marker', '.', 'MarkerSize', markerSize);

% Customize axes
ax = 1;
allH(ax).XLim            = [settings.nObservations(1)-1000, settings.nObservations(end)+500];
allH(ax).YLim            = [0 4.2];
allH(ax).XTick           = settings.nObservations;
allH(ax).XTickLabel      = settings.nObservations;
allH(ax).XAxis.FontSize  = 8;
allH(ax).YAxis.FontSize  = 8;
allH(ax).XAxis.TickLabelRotation = 90;
allH(ax).YTick           = 1:1:4;
allH(ax).YTickLabel      = 1:1:4;
allH(ax).XLabel.String   = 'Number of observations';
allH(ax).YLabel.String   = 'Time elapsed (s)';
allH(ax).XLabel.FontSize = 9;
allH(ax).YLabel.FontSize = 9;

allH(ax).XColor = [1 1 1];
allH(ax).YColor = [1 1 1];
allH(ax).Color  = [0 0 0];

print(fullfile(outDir, 'GraphicalAbstract_numObs.png'), '-dpng', '-r900');
close(fig);