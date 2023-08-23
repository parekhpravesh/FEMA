%% Plot results for simulation 001: effect of binning
% Load all results
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/fromTSD/';
load(fullfile(workDir, 'Results_experiment01.mat'));

%% Analyze bins
% For every bin, calculate the mean (over repeats) squared error of fixed
% effects from ground truth. Additionally, calculate the average 
% computational time for every bin, where the average is calculated across 
% repeats

% Initialize
allBins  = settings.allBins;
allMSE   = zeros(length(allBins), settings.nXvars, settings.nyVars);
avgTime  = zeros(length(allBins), 1);

% Calculate MSE for every bin for every y and every fixed effect parameter
for bins = 1:length(allBins)
    currResults  = results(:,bins);
    currResults  = vertcat(currResults{:});
    tmpEstimates = cat(3, currResults(:).beta_hat);

    % MSE is mean of the squared differences between estimates and
    % ground truth; the mean is taken after squaring the differences
    % across the repeats
    allMSE(bins, 1:settings.nXvars, :) = mean((tmpEstimates(2:end,:,:) - settings.GTruth.beta).^2, 3);

    % Average computational time for the bin
    avgTime(bins, 1) = mean(cat(1,currResults.elapsed));
end

%% Total sum of MSE across y variables
totalMSE     = sum(allMSE,3);
meanTotalMSE = mean(totalMSE,2);
minBinValue  = allBins(meanTotalMSE == min(meanTotalMSE));
MSE_20       = meanTotalMSE(allBins == 20);
MSE_min      = meanTotalMSE(allBins == minBinValue);

%% Computational time at 20 and min bin
cputime20  = avgTime(allBins == 20);
cputimeMin = avgTime(allBins == minBinValue);

%% Prepare plot
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 10]);
allH = tight_subplot(1, 2, 0.1, [0.14 0.02], [0.11, 0.02]);

% Plot total MSE
semilogx(allH(1), allBins, totalMSE, 'LineStyle', ':', 'LineWidth', 0.75, 'Color', [0 0 0], 'Marker', '.');
hold(allH(1), 'on');
allH(1).Box = 'off';

% Adjust X axis
allH(1).XLim = [1 max(allBins)];
allH(1).XAxis.TickValues = allBins;
allH(1).XAxis.MinorTick = 'off';
allH(1).XAxis.TickLabelRotation = 90;
allH(1).XAxis.FontSize = 8;
allH(1).XAxis.Label.String = 'Bins';
allH(1).XAxis.Label.FontSize = 12;
allH(1).XAxis.FontName = 'CMUBright';

% Adjust Y axis
allH(1).YAxis.FontSize = 10;
allH(1).YAxis.Label.String = 'Total MSE across imaging variables';
allH(1).YAxis.Label.FontSize = 12;
allH(1).YAxis.TickLabelFormat = '%.3f';
allH(1).YAxis.FontName = 'CMUBright';

% Average across all five X variables
semilogx(allH(1), allBins, meanTotalMSE, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [217 95 2]./255, 'Marker', '.', 'MarkerSize', 14);

% Get min and max of y
Ylimits_H1 = allH(1).YLim;

% Vertical line at bin 20
plot(allH(1), repmat(20, 100, 1), linspace(Ylimits_H1(1), MSE_20, 100), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [117 112 179]./255);
plot(allH(1), linspace(1,20,100), repmat(MSE_20, 100, 1), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [117 112 179]./255);
allH(1).YLim = Ylimits_H1;

% Vertical line at min bin
plot(allH(1), repmat(minBinValue, 100, 1), linspace(Ylimits_H1(1), MSE_min, 100), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [27 158 119]./255);
plot(allH(1), linspace(1, minBinValue, 100), repmat(MSE_min, 100, 1), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [27 158 119]./255);

% Add text
text(allH(1), 20,          MSE_20+0.0005,  num2str(MSE_20, '%.5f'),  'FontSize', 7, 'HorizontalAlignment', 'center', 'FontName', 'CMUBright');
text(allH(1), minBinValue, MSE_min+0.0005, num2str(MSE_min, '%.5f'), 'FontSize', 7, 'HorizontalAlignment', 'center', 'FontName', 'CMUBright');

% Label figure as a
text(allH(1), 0.13, Ylimits_H1(2), 'a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% Plot computational time required for bins
semilogx(allH(2), allBins, avgTime, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', [217 95 2]./255, 'Marker', '.', 'MarkerSize', 14);
hold(allH(2), 'on');
allH(2).Box = 'off';

% Adjust X axis
allH(2).XLim = [1 max(allBins)];
allH(2).XAxis.TickValues = allBins;
allH(2).XAxis.MinorTick = 'off';
allH(2).XAxis.TickLabelRotation = 90;
allH(2).XAxis.FontSize = 8;
allH(2).XAxis.Label.String = 'Bins';
allH(2).XAxis.Label.FontSize = 12;
allH(2).XAxis.FontName = 'CMUBright';

% Adjust Y axis
allH(2).YAxis.FontSize = 10;
allH(2).YAxis.Label.String = 'Average CPU time (s)';
allH(2).YAxis.Label.FontSize = 12;
allH(2).YAxis.FontName = 'CMUBright';

% Vertical line at bin 20
Ylimits_H2 = allH(2).YLim;
plot(allH(2), repmat(20, 100, 1), linspace(Ylimits_H2(1), cputime20, 100), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [117 112 179]./255);
plot(allH(2), linspace(1,20,100), repmat(cputime20, 100, 1), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [117 112 179]./255);
allH(2).YLim = Ylimits_H2;

% Vertical line at min bin
plot(allH(2), repmat(minBinValue, 100, 1),   linspace(Ylimits_H2(1), cputimeMin, 100), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [27 158 119]./255);
plot(allH(2), linspace(1, minBinValue, 100), repmat(cputimeMin, 100, 1), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', [27 158 119]./255);

% Add text
text(allH(2), 20,          cputime20+2,  num2str(cputime20, '%.2f'),  'FontSize', 7, 'HorizontalAlignment', 'right', 'FontName', 'CMUBright');
text(allH(2), minBinValue, cputimeMin+2, num2str(cputimeMin,'%.2f'),  'FontSize', 7, 'HorizontalAlignment', 'right', 'FontName', 'CMUBright');

% Label figure as b
text(allH(2), 0.22, Ylimits_H2(2), 'b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% Save this figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/Results01-EffectBins.png', '-dpng', '-r900');
close(fig);