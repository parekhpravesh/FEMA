%% Summarize and plot results from simulation 001 - parameter recovery at bin size 20
%% Load all results
workDir = '/Applications/Projects/2023-01-31_FEMAExperiments/fromTSD/';
outDir  = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-08-11_FiguresRedone/simulation-GTruth';
load(fullfile(workDir, 'Results_experiment01.mat'));

%% Analyze bin == 20
bin     = 20;
binLoc  = find(settings.allBins == bin);

currResults  = results(:,binLoc);
currResults  = vertcat(currResults{:});

overallEstimates_FFX      = cat(3, currResults.beta_hat);
overallEstimates_RFX      = cat(3, currResults.sig2mat);
overallEstimates_sig2tvec = cat(3, currResults.sig2tvec);

toUse_estimates_RFX = overallEstimates_RFX .* overallEstimates_sig2tvec;

%% Fixed effects
for ffx = 1:5
    fig  = figure('Units', 'centimeters', 'Position', [10 10 5 5]);
    allH = tight_subplot(1, 1, 0, [0.22 0.03], [0.25 0.038]);
    hold(allH(1), 'on');
    scatter(allH(1), settings.GTruth.beta(ffx,:), squeeze(overallEstimates_FFX(ffx+1,:,:)),         'Marker', '.', 'MarkerEdgeColor', [231,212,232]./255, 'SizeData', 5);
    scatter(allH(1), settings.GTruth.beta(ffx,:), mean(squeeze(overallEstimates_FFX(ffx+1,:,:)),2), 'Marker', '.', 'MarkerFaceColor', [118,42,131]./255, 'MarkerEdgeColor', [118,42,131]./255, 'SizeData', 5);
    allH(1).XLim              = [-0.2 0.2];
    allH(1).YLim              = [-0.2 0.2];
    allH(1).XTick             = -0.2:0.1:0.2;
    allH(1).YTick             = -0.2:0.1:0.2;
    allH(1).XTickLabelMode    = 'auto';
    allH(1).YTickLabelMode    = 'auto';
    allH(1).XAxis.FontSize    = 8;
    allH(1).YAxis.FontSize    = 8;
    allH(1).YLabel.String     = 'Estimated \beta';
    allH(1).XLabel.String     = 'Simulated \beta';
    allH(1).XAxis.FontName    = 'CMUBright';
    allH(1).YAxis.FontName    = 'CMUBright';
    allH(1).XLabel.FontSize   = 11;
    allH(1).YLabel.FontSize   = 11;
    allH(1).XAxis.TickLabelRotation = 0;
    print(fullfile(outDir, ['FFX-', num2str(ffx, '%02d'), '.png']), '-dpng', '-r900');
    close(fig);
end

%% Random effects
for rfx = 1:3
    fig  = figure('Units', 'centimeters', 'Position', [10 10 5 5]);
    allH = tight_subplot(1, 1, 0, [0.22 0.03], [0.25 0.038]);
    hold(allH(1), 'on');
    scatter(allH(1), settings.GTruth.sig2mat_true(rfx,:), squeeze(toUse_estimates_RFX(rfx,:,:)),         'Marker', '.', 'MarkerEdgeColor', [217,240,211]./255, 'SizeData', 5);
    scatter(allH(1), settings.GTruth.sig2mat_true(rfx,:), mean(squeeze(toUse_estimates_RFX(rfx,:,:)),2), 'Marker', '.', 'MarkerFaceColor', [27,120,55]./255, 'MarkerEdgeColor', [27,120,55]./255, 'SizeData', 5);
    if rfx == 3
        allH(1).XLim          = [0.2 0.8];
        allH(1).YLim          = [0.2 0.8];
    else
        allH(1).XLim          = [0.0 0.8];
        allH(1).YLim          = [0.0 0.8];
    end
    allH(1).XTick             = 0.0:0.2:0.8;
    allH(1).YTick             = 0.0:0.2:0.8;
    allH(1).XTickLabelMode    = 'auto';
    allH(1).YTickLabelMode    = 'auto';
    allH(1).XAxis.FontSize    = 8;
    allH(1).YAxis.FontSize    = 8;
    allH(1).YLabel.String     = 'Estimated \sigma^2';
    allH(1).XLabel.String     = 'Simulated \sigma^2';
    allH(1).XAxis.FontName    = 'CMUBright';
    allH(1).YAxis.FontName    = 'CMUBright';
    allH(1).XLabel.FontSize   = 11;
    allH(1).YLabel.FontSize   = 11;
    allH(1).XAxis.TickLabelRotation = 0;
    print(fullfile(outDir, ['RFX-', num2str(rfx, '%02d'), '.png']), '-dpng', '-r900');
    close(fig);
end