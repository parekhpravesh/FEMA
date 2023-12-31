%% Plot results for simulation 004: type I error rate
%% Load results
rootDir = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/export_TSD/';
toLoad  = 'Results_experiment04.mat';
load(fullfile(rootDir, toLoad));

%% Locations for different variables
loc_FEMA_FP   = find(~(cellfun(@isempty, regexpi(results.Properties.VariableNames, '^FP_FEMA_y'))));
loc_fitlme_FP = find(~(cellfun(@isempty, regexpi(results.Properties.VariableNames, '^FP_fitlmematrix_y'))));

%% Colour scheme
col_FEMA        = [27, 158, 119]./255;
col_fitlmatrix  = [217, 95,   2]./255;

%% Plot results
fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
allH = tight_subplot(3, 1, [0.08 0.01], [0.215 0.01], [0.07 0.01]);

hold(allH(1), 'on');
hold(allH(2), 'on');

% Customize axis size
% Move severally right
allH(3).Position(1) = allH(3).Position(1) + 0.7;

% Shorten
allH(3).Position(3) = allH(3).Position(3) - 0.7;

% Align top
allH(3).Position(2) = allH(2).Position(2);

% Shrink axis 2
allH(2).Position(3) = allH(2).Position(3) - 0.35;

% Overall need to expand all axes vertically
allH(1).Position(4) = allH(1).Position(4) + 0.18;
allH(2).Position(4) = allH(2).Position(4) + 0.18;
allH(3).Position(4) = allH(3).Position(4) + 0.18;

% Now move their base
allH(1).Position(2) = allH(1).Position(2) - 0.18;
allH(2).Position(2) = allH(2).Position(2) - 0.4;
allH(3).Position(2) = allH(3).Position(2) - 0.4;

% False positives for each y variable, averaged across repeats
means_fitlme = mean(results{:, loc_fitlme_FP});
means_FEMA   = mean(results{:, loc_FEMA_FP});

% Overall maximum
maxAll1 = ceil(max([results{:, loc_fitlme_FP}; results{:, loc_FEMA_FP}], [], 'all'));

% Dashed line at 5
l3 = plot(allH(1), 0:31, repmat(5, 32, 1), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
l4 = plot(allH(2), 0:21, repmat(5, 22, 1), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);

for bars = 1:30
    % Add scatter plots within each bar
    d1 = scatter(allH(1), repmat(bars-0.15, 100, 1), results{:,loc_fitlme_FP(bars)} + rand(100,1)./2, 'o', 'MarkerEdgeColor', col_fitlmatrix, 'SizeData', 25, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.65, 'MarkerFaceAlpha', 0.75);
    d2 = scatter(allH(1), repmat(bars+0.15, 100, 1), results{:,loc_FEMA_FP(bars)}   + rand(100,1)./2, 'o', 'MarkerEdgeColor', col_FEMA,       'SizeData', 25, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.65, 'MarkerFaceAlpha', 0.75);
    d1.SizeData = d1.SizeData/8;
    d2.SizeData = d2.SizeData/8;

    bb1           = bar(allH(1), bars-0.15, ceil(means_fitlme(bars)), 0.2, 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 0.75, 'BarLayout','grouped');
    bb2           = bar(allH(1), bars+0.15, ceil(means_FEMA(bars)),   0.2, 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 0.75, 'BarLayout','grouped');
    bb1.CData     = col_fitlmatrix;
    bb2.CData     = col_FEMA;
    bb1.EdgeColor = col_fitlmatrix;
    bb2.EdgeColor = col_FEMA;
end

count = 1;
for bars = 31:50
    % Add scatter plots within each bar
    d1 = scatter(allH(2), repmat(count-0.15, 100, 1), results{:,loc_fitlme_FP(bars)} + rand(100,1)./2, 'o', 'MarkerEdgeColor', col_fitlmatrix, 'SizeData', 25, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.65, 'MarkerFaceAlpha', 0.75);
    d2 = scatter(allH(2), repmat(count+0.15, 100, 1), results{:,loc_FEMA_FP(bars)}   + rand(100,1)./2, 'o', 'MarkerEdgeColor', col_FEMA,       'SizeData', 25, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.65, 'MarkerFaceAlpha', 0.75);
    d1.SizeData = d1.SizeData/8;
    d2.SizeData = d2.SizeData/8;

    bb3           = bar(allH(2), count-0.15, ceil(means_fitlme(bars)), 0.2, 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 0.75, 'BarLayout','grouped');
    bb4           = bar(allH(2), count+0.15, ceil(means_FEMA(bars)),   0.2, 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 0.75, 'BarLayout','grouped');
    bb3.CData     = col_fitlmatrix;
    bb4.CData     = col_FEMA;
    bb3.EdgeColor = col_fitlmatrix;
    bb4.EdgeColor = col_FEMA;

    count = count + 1;
end

% Customize axes
allH(1).XLim             = [0.5 30.5];
allH(1).YLim             = [0 maxAll1+1];
allH(1).XTick            = 1:1:30;
allH(1).XTickLabel       = 1:1:30;
allH(1).XAxis.FontSize   = 8;
allH(1).XLabel.FontSize  = 10;
allH(1).YTick            = 5:5:maxAll1+1;
allH(1).YTickLabel       = 5:5:maxAll1+1;
allH(1).YGrid            = 'on';
allH(1).YMinorGrid       = 'on';
allH(1).XAxis.TickLength = [0 0];
allH(1).YLabel.String    = 'Number of FPs (out of 100)';
allH(1).YLabel.FontSize  = 10;
allH(1).XLabel.String    = 'Imaging variable';
allH(1).XAxis.FontName   = 'CMUBright';
allH(1).YAxis.FontName   = 'CMUBright';
allH(1).XAxis.TickLabelRotation = 0;

allH(2).XLim             = [0.5 20.5];
allH(2).YLim             = [0 maxAll1+1];
allH(2).XTick            = 1:1:20;
allH(2).XTickLabel       = 31:1:50;
allH(2).XAxis.FontSize   = 8;
allH(2).XLabel.FontSize  = 10;
allH(2).YTick            = 5:5:maxAll1+1;
allH(2).YTickLabel       = 5:5:maxAll1+1;
allH(2).YGrid            = 'on';
allH(2).YMinorGrid       = 'on';
allH(2).XAxis.TickLength = [0 0];
allH(2).YLabel.String    = 'Number of FPs (out of 100)';
allH(2).YLabel.FontSize  = 10;
allH(2).XLabel.String    = 'Imaging variable';
allH(2).XAxis.FontName   = 'CMUBright';
allH(2).YAxis.FontName   = 'CMUBright';
allH(2).XAxis.TickLabelRotation = 0;

% Turn on legend
ll1 = legend(allH(1), [bb2, bb1, l3], {'FEMA', 'fitlmematrix', '5% false positives (FPs)'}, 'Box', 'off', 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontName', 'CMUBright');
ll2 = legend(allH(2), [bb4, bb3, l4], {'FEMA', 'fitlmematrix', '5% false positives (FPs)'}, 'Box', 'off', 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontName', 'CMUBright');

% Touch the legends just a bit so that the spacing goes back to what it was
ll1.Position(1) = ll1.Position(1) + eps;
ll2.Position(1) = ll2.Position(1) + eps;

% Move the first legend down
ll1.Position(2) = ll1.Position(2) - 0.1;

% Do the third plot - overall total number of false positives, averaged across repetitions
hold(allH(3), 'on');

bb = bar(allH(3), 1, [ceil(mean(results.FP_fitlmematrix)); ceil(mean(results.FP_FEMA))], 'grouped', 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 1);
bb(1).CData     = col_fitlmatrix;
bb(2).CData     = col_FEMA;
bb(1).EdgeColor = col_fitlmatrix;
bb(2).EdgeColor = col_FEMA;

% Add scatter plots within each
scatter(allH(3), repmat(1-0.15, 100, 1), results.FP_fitlmematrix, 'o', 'MarkerEdgeColor', col_fitlmatrix, 'SizeData', 15, 'jitter', 'on', 'jitterAmount', 0.025);
scatter(allH(3), repmat(1+0.15, 100, 1), results.FP_FEMA,         'o', 'MarkerEdgeColor', col_FEMA,       'SizeData', 15, 'jitter', 'on', 'jitterAmount', 0.025);

% Dashed line at 250
l3 = plot(allH(3), 0:2, repmat(250, 3, 1), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);

% Customize axes
maxAll2                 = ceil(max([results.FP_fitlmematrix; results.FP_FEMA], [], 'all'));
allH(3).XLim            = [0.5 1.5];
allH(3).YLim            = [0 maxAll2+9];
allH(3).XTick           = [];
allH(3).XTickLabel      = [];
allH(3).YTick           = 0:50:maxAll2+9;
allH(3).YTickLabel      = 0:50:maxAll2+9;
allH(3).YLabel.String   = 'Number of FPs across all imaging variables';
allH(3).YLabel.FontSize = 8;
allH(3).XAxis.FontName  = 'CMUBright';
allH(3).YAxis.FontName  = 'CMUBright';

% Move legend one to be somewhere in the centre
ll2.Position(2) = ll2.Position(2) - 0.1;
ll2.Position(1) = ll1.Position(1);
% ll2.Position(1) = ll2.Position(1) + 0.15;

% Align bottom
pause(1);
allH(3).Position(2) = allH(2).Position(2);

% Put text a and b
t1 = text(allH(1), -1.7,  maxAll1+1, 'a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');
t2 = text(allH(2), -1.8,  maxAll1+1, 'b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');
t3 = text(allH(3), -0.01, maxAll2+7, 'c)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');

% Save the figure
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/Results04-TypeI.png', '-dpng', '-r900');
close(fig);

%% Deprecated code for 20 imaging variables
% %% Plot results
% fig  = figure('Units', 'centimeters', 'Position', [10 10 16 10]);
% allH = tight_subplot(1, 2, 0.2, [0.2 0.015], [0.07 0.01]);
% hold(allH(1), 'on');
% 
% % Customize axis size
% % Move severally right
% allH(2).Position(1) = allH(2).Position(1) + 0.2;
% 
% % Shorten
% allH(2).Position(3) = allH(2).Position(3) - 0.2;
% 
% % Align top
% allH(2).Position(4) = allH(1).Position(4);
% 
% % Expand axis 1
% allH(1).Position(3) = allH(1).Position(3) + 0.3;
% 
% % False positives for each y variable, averaged across repeats
% means_fitlme = mean(results{:, loc_fitlme_FP});
% means_FEMA   = mean(results{:, loc_FEMA_FP});
% 
% % Overall maximum
% maxAll1 = ceil(max([results{:, loc_fitlme_FP}; results{:, loc_FEMA_FP}], [], 'all'));
% 
% % Dashed line at 5
% l3 = plot(allH(1), 0:21, repmat(5, 22, 1), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
% 
% for bars = 1:20
% 
%     % Add scatter plots within each bar
%     d1 = scatter(allH(1), repmat(bars-0.15, 100, 1), results{:,loc_fitlme_FP(bars)} + rand(100,1)./2, 'o', 'MarkerEdgeColor', col_fitlmatrix, 'SizeData', 25, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.65, 'MarkerFaceAlpha', 0.75);
%     d2 = scatter(allH(1), repmat(bars+0.15, 100, 1), results{:,loc_FEMA_FP(bars)}   + rand(100,1)./2, 'o', 'MarkerEdgeColor', col_FEMA,       'SizeData', 25, 'jitter', 'on', 'jitterAmount', 0.1, 'LineWidth', 0.5, 'MarkerEdgeAlpha', 0.65, 'MarkerFaceAlpha', 0.75);
%     d1.SizeData = d1.SizeData/8;
%     d2.SizeData = d2.SizeData/8;
% 
%     bb1           = bar(allH(1), bars-0.15, ceil(means_fitlme(bars)), 0.2, 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 0.75, 'BarLayout','grouped');
%     bb2           = bar(allH(1), bars+0.15, ceil(means_FEMA(bars)),   0.2, 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 0.75, 'BarLayout','grouped');
%     bb1.CData     = col_fitlmatrix;
%     bb2.CData     = col_FEMA;
%     bb1.EdgeColor = col_fitlmatrix;
%     bb2.EdgeColor = col_FEMA;
% 
% end
% 
% % Customize axes
% allH(1).XLim             = [0.5 20.5];
% allH(1).YLim             = [0 maxAll1+0.5];
% allH(1).XTick            = 1:1:20;
% allH(1).XTickLabel       = 1:1:20;
% allH(1).XAxis.FontSize   = 8;
% allH(1).XLabel.FontSize  = 10;
% allH(1).YTick            = 5:5:maxAll1;
% allH(1).YTickLabel       = 5:5:maxAll1;
% allH(1).YGrid            = 'on';
% allH(1).YMinorGrid       = 'on';
% allH(1).XAxis.TickLength = [0 0];
% allH(1).YLabel.String    = 'Number of FPs (out of 100)';
% allH(1).YLabel.FontSize  = 10;
% allH(1).XLabel.String    = 'Imaging variable';
% allH(1).XAxis.FontName   = 'CMUBright';
% allH(1).YAxis.FontName   = 'CMUBright';
% allH(1).XAxis.TickLabelRotation = 0;
% 
% % Turn on legend
% ll = legend(allH(1), [bb2, bb1, l3], {'FEMA', 'fitlmematrix', '5% false positives (FPs)'}, 'Box', 'off', 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontName', 'CMUBright');
% 
% % Do the second plot - overall total number of false positives, averaged across repetitions
% hold(allH(2), 'on');
% 
% bb = bar(allH(2), 1, [ceil(mean(results.FP_fitlmematrix)); ceil(mean(results.FP_FEMA))], 'grouped', 'Horizontal', 'off', 'FaceColor', 'none', 'LineWidth', 1);
% bb(1).CData     = col_fitlmatrix;
% bb(2).CData     = col_FEMA;
% bb(1).EdgeColor = col_fitlmatrix;
% bb(2).EdgeColor = col_FEMA;
% 
% % Add scatter plots within each
% scatter(allH(2), repmat(1-0.15, 100, 1), results.FP_fitlmematrix, 'o', 'MarkerEdgeColor', col_fitlmatrix, 'SizeData', 15, 'jitter', 'on', 'jitterAmount', 0.025);
% scatter(allH(2), repmat(1+0.15, 100, 1), results.FP_FEMA,         'o', 'MarkerEdgeColor', col_FEMA,       'SizeData', 15, 'jitter', 'on', 'jitterAmount', 0.025);
% 
% % Dashed line at 100
% l3 = plot(allH(2), 0:2, repmat(100, 3, 1), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
% 
% % Customize axes
% maxAll2                 = ceil(max([results.FP_fitlmematrix; results.FP_FEMA], [], 'all'));
% allH(2).XLim            = [0.5 1.5];
% allH(2).YLim            = [0 maxAll2+9];
% allH(2).XTick           = [];
% allH(2).XTickLabel      = [];
% allH(2).YTick           = 0:20:maxAll2+9;
% allH(2).YTickLabel      = 0:20:maxAll2+9;
% allH(2).YLabel.String   = 'Number of FPs across all imaging variables';
% allH(2).YLabel.FontSize = 10;
% allH(2).XAxis.FontName  = 'CMUBright';
% allH(2).YAxis.FontName  = 'CMUBright';
% 
% % Move legend one to be somewhere in the centre
% ll.Position(2) = ll.Position(2) - 0.18;
% ll.Position(1) = ll.Position(1) + 0.125;
% 
% % Align bottom
% pause(1);
% allH(2).Position(2) = allH(1).Position(2);
% 
% % Put text a and b
% t1 = text(allH(1), -1.5,  maxAll1+0.5, 'a)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');
% t2 = text(allH(2), -0.01, maxAll2+9,   'b)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'CMUBright');