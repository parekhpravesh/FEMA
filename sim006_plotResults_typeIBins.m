%% Plot results for simulation 006: type I error rate as a function of bins
%% Load results
rootDir = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/export_TSD/';
toLoad  = 'Results_experiment06.mat';
load(fullfile(rootDir, toLoad));
allBins = 1:30;

%% Colour scheme
col_FEMA = [27, 158, 119]./255;

%% Compile average false positives in each bin
allAverageFP = zeros(length(allBins), 1);
for bins = 1:length(allBins)
    allAverageFP(bins,1) = ceil(mean(results.FP_FEMA(results.BinSize == allBins(bins))));
end

%% Plot results
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 9]);
allH = tight_subplot(1, 1, [0.08 0.01], [0.1 0.01], [0.1 0.01]);

hold(allH(1), 'on');

% Overall maximum
maxAll1 = ceil(max(results.FP_FEMA));

% Dashed line at 2500
l3 = plot(allH(1), 0:31, repmat(2500, 32, 1), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);

for bins = 1:30
    % Which rows?
    toRows = results.BinSize == allBins(bins);

    % Add scatter plots within each bar
    d2 = scatter(allH(1), repmat(bins, 100, 1), results.FP_FEMA(toRows) + rand(100,1)./2, 'o', 'MarkerEdgeColor', col_FEMA, 'SizeData', 25, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.65, 'MarkerFaceAlpha', 0.75);
    d2.SizeData = d2.SizeData/8;

    b2 = bar(allH(1), bins, ceil(mean(results.FP_FEMA(toRows))), 0.5, 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 0.75, 'BarLayout','grouped');
    % b2.CData     = col_FEMA;
    % b2.EdgeColor = col_FEMA;
end

% Overall line for average false positives
% plot(allH(1), allBins, allAverageFP, 'LineStyle', '-', 'LineWidth', 1);

% Customize axes
allH(1).XLim             = [0.5 30.5];
allH(1).YLim             = [0 maxAll1+50];
allH(1).XTick            = 1:1:30;
allH(1).XTickLabel       = allBins;
allH(1).XAxis.FontSize   = 8;
allH(1).XLabel.FontSize  = 10;
allH(1).YTick            = 500:500:maxAll1+1;
allH(1).YTickLabel       = 500:500:maxAll1+1;
allH(1).XAxis.TickLength = [0 0];
allH(1).YLabel.String    = 'Number of FPs across imaging variables';
allH(1).YLabel.FontSize  = 10;
allH(1).XLabel.String    = 'Bin Value';
allH(1).XAxis.FontName   = 'CMUBright';
allH(1).YAxis.FontName   = 'CMUBright';
allH(1).XAxis.TickLabelRotation = 0;

%% Save the figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results06-TypeIBins.png', '-dpng', '-r900');
close(fig);