%% Plot results of simulation 002: comparison with fitlme
%% Load all results
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/fromTSD/';
load(fullfile(workDir, 'Results_experiment02.mat'));

%% Sample size of 10000
sampSize = 10000;
sampLoc  = find(settings.nObservations == sampSize);

%% Subset results
subsetResults_FEMA  = results_FEMA{sampLoc};
subsetResults_LME   = results_fitlme{sampLoc};

%% Plot all parameter estimates with SE
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 18], 'PaperUnits','centimeters', 'PaperSize',[16 16]);
allH = tight_subplot(1, 5, [0.0 0.01], [0.12 0.03], [0.07 0.01]);

% Color scheme
col_fitlme  = [217,95,2]./255;
col_FEMA    = [27,158,119]./255;

% Locations to plot at
xx1 = 0.7:1:49.7;
xx2 = 1.3:1:50.3;

% Hold everything
hold(allH(:), 'on');

% Settings
markerType  = 'o';
markerSize  = 3;
LineStyle   = 'none';
LineWidth   = 0.5;
capSize     = 4;

% Plot
for ax = 1:5

    % Actual plotting
    errorbar(allH(ax), subsetResults_FEMA.beta_hat_FEMA(ax,  1:end), xx2, subsetResults_FEMA.beta_se_FEMA(ax,  1:end), 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    errorbar(allH(ax), subsetResults_LME.beta_hat_fitlme(ax, 1:end), xx1, subsetResults_LME.beta_se_fitlme(ax, 1:end), 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);

    % Apply axes settings
    allH(ax).YLim            = [0 51];
    allH(ax).Box             = 'off';
    allH(ax).XLim            = [-0.25 0.25];
    allH(ax).YTickLabel      = [];
    allH(ax).YMinorGrid      = 'on';
    allH(ax).XLabel.String   = 'Estimated \beta';
    allH(ax).Title.String    = ['X', num2str(ax)];
    allH(ax).Title.FontWeight = 'normal';
    allH(ax).Title.FontName  = 'CMUBright';
    allH(ax).XAxis.FontSize  = 8;
    allH(ax).XLabel.FontSize = 11;
    allH(ax).YLabel.FontSize = 11;
    allH(ax).Title.FontSize  = 11;
    allH(ax).XTick           = -0.2:0.2:0.2;
    allH(ax).XTickLabel      = -0.2:0.2:0.2;
    allH(ax).XAxis.FontName  = 'CMUBright';
    allH(ax).YAxis.FontName  = 'CMUBright';
end

allH(1).YLabel.String   = 'Imaging variables';
allH(1).YTick           = 5:5:50;
allH(1).YTickLabel      = num2str((5:5:50)', '%02d');

% Overall legend
ll = legend(allH(2), {'FEMA', 'fitlme'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off', 'FontName', 'CMUBright');
ll.Position(2) = ll.Position(2) - 0.11;
ll.Position(1) = allH(3).Position(1) - 0.1;

% Save this figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/Results03-Comparison_MATLAB_FEMA_Fixed.png', '-dpng', '-r900');
close(fig);

%% Comparison of random effects
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 18], 'PaperUnits','centimeters', 'PaperSize',[16 16]);
allH = tight_subplot(1, 3, [0.0 0.01], [0.12 0.03], [0.07 0.01]);

% Color scheme
col_fitlme  = [217,95,2]./255;
col_FEMA    = [27,158,119]./255;

% Locations to plot at
xx1 = 0.7:1:49.7;
xx2 = 1.3:1:50.3;

% Hold all plots
hold(allH(:), 'on');

% Plot
for ax = 1:3
    errorbar(allH(ax), subsetResults_FEMA.sig2mat_FEMA(ax,   1:end),    xx2, subsetResults_FEMA.sig2Low_FEMA(ax,  1:end),    subsetResults_FEMA.sig2Upp_FEMA(ax,  1:end),    'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
    errorbar(allH(ax), subsetResults_LME.sig2mat_fitlme(ax,  1:end).^2, xx1, subsetResults_LME.sig2Low_fitlme(ax, 1:end).^2, subsetResults_LME.sig2Upp_fitlme(ax, 1:end).^2, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);

    % Apply axes settings
    allH(ax).YLim               = [0 51];
    allH(ax).Box                = 'off';
    allH(ax).XLim               = [-0.1 1.8];
    allH(ax).YTickLabel         = [];
    allH(ax).YMinorGrid         = 'on';
    allH(ax).XTick              = 0:0.4:1.8;
    allH(ax).XTickLabel         = 0:0.4:1.8;
    allH(ax).YAxis.FontSize     = 8;
    allH(ax).XAxis.FontSize     = 8;
    allH(ax).YLabel.FontSize    = 11;
    allH(ax).XLabel.FontSize    = 11;
    allH(ax).Title.FontSize     = 11;
    allH(ax).Title.FontWeight   = 'normal';
    allH(ax).XLabel.String      = 'Estimated \sigma^2';
    allH(ax).XAxis.FontName     = 'CMUBright';
    allH(ax).YAxis.FontName     = 'CMUBright';
end

allH(1).YLabel.String   = 'Imaging variables';
allH(1).Title.String    = 'Family effect (\itF\rm)';
allH(2).Title.String    = 'Subject effect (\itS\rm)';
allH(3).Title.String    = 'Unmodeled variance (\itE\rm)';

allH(1).Title.FontName  = 'CMUBright';
allH(2).Title.FontName  = 'CMUBright';
allH(3).Title.FontName  = 'CMUBright';

allH(1).YTick           = 5:5:50;
allH(1).YTickLabel      = num2str((5:5:50)', '%02d');

% Overall legend
ll              = legend(allH(2), {'FEMA', 'fitlme'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off', 'FontName', 'CMUBright');
ll.Position(2)  = ll.Position(2) - 0.11;

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/Results03-Comparison_MATLAB_FEMA_RFX.png', '-dpng', '-r900');
close(fig);