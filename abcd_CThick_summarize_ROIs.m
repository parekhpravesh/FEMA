%% Plot results for DK40 - MATLAB vs. FEMA
% Get the correct version of tight_subplot (which returns axes handles)
addpath('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/code/cmig_tools_internal-beta/showSurf');

% Get data
dirResults = '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/CorticalThickness_DK_FSE_nn';
load(fullfile(dirResults, 'Results.mat'), 'beta_*', 'covariates', 'covarNames', 'logpmat', 'lowVar_MATLAB', 'RandomEffects', 'roiNames', 'sig2*', 'UppVar_MATLAB', 'varComp_MATLAB', 'zmat', 'mdl');

% Identify column of interest - should be third column: intercept, age, age
% delta
locAgeDelta = strcmpi(covarNames, 'ageDelta');
locAge      = strcmpi(covarNames, 'ageBaseline');

%% Compile p values 
pValues_MATLAB = zeros(length(covarNames), 68);
for rois = 1:68
    pValues_MATLAB(:,rois) = mdl{rois}.Coefficients.pValue;
end

%% Return to actual p values - FEMA
pValues_FEMA = normcdf(-abs(zmat))*2;

%% New roi names
roiNames_Use = strrep(roiNames, 'smri_thick_cdk_', '');

%% Identify left and right hemispheres
locLeft     = find(~cellfun(@isempty, regexpi(roiNames_Use, 'lh$')));
locRight    = find(~cellfun(@isempty, regexpi(roiNames_Use, 'rh$')));

%% Make a summary table of ROIs, beta values, and random effects
results = cell2table([roiNames_Use', num2cell(...
[beta_MATLAB(locAgeDelta, :)', beta_se_MATLAB(locAgeDelta, :)', pValues_MATLAB(locAgeDelta, :)',   ...
 beta_hat(locAgeDelta, :)',    beta_se(locAgeDelta, :)',        pValues_FEMA(locAgeDelta, :)',     ...
 beta_MATLAB(locAge, :)',      beta_se_MATLAB(locAge, :)',      pValues_MATLAB(locAge, :)',        ...
 beta_hat(locAge, :)',         beta_se(locAge, :)',             pValues_FEMA(locAge, :)',          ...
 varComp_MATLAB(1,:)',         lowVar_MATLAB(1,:)',             UppVar_MATLAB(1,:)',               ...
 (sig2mat(1, :).*sig2tvec)',   sig2Low(1,:)',                   sig2Upp(1, :)',                    ...
 varComp_MATLAB(2,:)',         lowVar_MATLAB(2,:)',             UppVar_MATLAB(2,:)',               ...
 (sig2mat(2, :).*sig2tvec)',   sig2Low(2,:)',                   sig2Upp(2, :)',                    ...
 varComp_MATLAB(3,:)',         lowVar_MATLAB(3,:)',             UppVar_MATLAB(3,:)',               ...
 (sig2mat(3, :).*sig2tvec)',   sig2Low(3,:)',                   sig2Upp(3, :)'])],                 ...
 'VariableNames', {'ROIName', 'Beta_ageDelta_MATLAB', 'SE_ageDelta_MATLAB', 'pValue_ageDelta_MATLAB', 'Beta_ageDelta_FEMA', 'SE_ageDelta_FEMA', 'pValue_ageDelta_FEMA', ...
                   'Beta_age_MATLAB', 'SE_age_MATLAB', 'pValue_age_MATLAB', 'Beta_age_FEMA', 'SE_age_FEMA', 'pValue_age_FEMA', ...
                   'varFamily_MATLAB', 'lowFamily_MATLAB', 'uppFamily_MATLAB', 'varFamily_FEMA', 'lowFamily_FEMA', 'uppFamily_FEMA', ...
                   'varSubject_MATLAB', 'lowSubject_MATLAB', 'uppSubject_MATLAB', 'varSubject_FEMA', 'lowSubject_FEMA', 'uppSubject_FEMA', ... 
                   'varError_MATLAB', 'lowError_MATLAB', 'uppError_MATLAB', 'varError_FEMA', 'lowError_FEMA', 'uppError_FEMA'});
writetable(results, '/home/pparekh/analyses/2023-02-17_FEMA-ABCD/CorticalThickness_DK_FSE_nn/Summary_DK40.xlsx');


%% Deprecated code for plotting that was previously used (prior to the summary section)
% %% Plot all parameter estimates
% fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16], 'PaperUnits','centimeters', 'PaperSize', [16 16]);
% allH = tight_subplot(1, 2, [0.0 0.04], [0.12 0.03], [0.12 0.02]);
% 
% % Color scheme
% col_fitlme  = [217,95,2]./255;
% col_FEMA    = [27,158,119]./255;
% 
% % Locations to plot at
% xx1 = 0.7:1:34.7;
% xx2 = 1.3:1:35.3;
% 
% % Hold everything
% hold(allH(:), 'on');
% 
% % Settings
% markerType  = 'o';
% markerSize  = 3;
% LineStyle   = 'none';
% LineWidth   = 0.5;
% capSize     = 4;
% 
% alpha = 0.05/68;
% 
% % Plot
% 
% % Slowly loop over all LH ROIs
% for rois = locLeft
%     % Determine if statistically significant or not - MATLAB
%     if pValues_MATLAB(locAgeDelta, locLeft(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(1), beta_MATLAB(locAgeDelta, locLeft(rois)), xx1(rois), beta_se_MATLAB(locAgeDelta, locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(1), beta_MATLAB(locAgeDelta, locLeft(rois)), xx1(rois), beta_se_MATLAB(locAgeDelta, locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% 
%     % Determine if statistically significant or not - FEMA
%     if pValues_FEMA(locAgeDelta, locLeft(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(1), beta_hat(locAgeDelta,    rois), xx2(rois), beta_se(locAgeDelta,        rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(1), beta_hat(locAgeDelta,    rois), xx2(rois), beta_se(locAgeDelta,        rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% end
% 
% % Slowly loop over all RH rois
% for rois = 1:length(locRight)
%     % Determine if statistically significant or not - MATLAB
%     if pValues_MATLAB(locAgeDelta, locRight(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(2), beta_MATLAB(locAgeDelta, locRight(rois)), xx1(rois), beta_se_MATLAB(locAgeDelta, locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(2), beta_MATLAB(locAgeDelta, locRight(rois)), xx1(rois), beta_se_MATLAB(locAgeDelta, locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% 
%     % Determine if statistically significant or not - FEMA
%     if pValues_FEMA(locAgeDelta, locRight(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(2), beta_hat(locAgeDelta,    locRight(rois)), xx2(rois), beta_se(locAgeDelta,        locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(2), beta_hat(locAgeDelta,    locRight(rois)), xx2(rois), beta_se(locAgeDelta,        locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% end
% 
% for ax = 1:2
%     % Apply axes settings
%     allH(ax).YLim            = [0 35];
%     allH(ax).Box             = 'off';
%     allH(ax).XLim            = [-0.026 0.026];
%     allH(ax).YTick           = 1:34;
%     allH(ax).YGrid           = 'on';
%     % allH(ax).YTickLabel      = [];
% %     allH(ax).YMinorGrid      = 'on';
%     % allH(ax).XMinorGrid      = 'off';
%     allH(ax).XLabel.String   = 'Estimated \beta age_{delta}';
%     % allH(ax).Title.String    = ['X', num2str(ax)];
%     allH(ax).Title.FontWeight = 'normal';
%     allH(ax).XAxis.FontSize  = 8;
%     allH(ax).XLabel.FontSize = 11;
%     allH(ax).YLabel.FontSize = 11;
%     allH(ax).Title.FontSize  = 11;
%     allH(ax).XTick           = -0.025:0.01:0.025;
%     allH(ax).XTickLabel      = -0.025:0.01:0.025;
% end
% 
% % allH(1).YLabel.String   = 'Imaging variables';
% allH(1).YTick           = 1:34;
% allH(1).YTickLabel      = strrep(roiNames_Use(locLeft), 'lh', '');
% % allH(1).TickLabelInterpreter = 'latex';
% 
% allH(1).Title.String    = 'Left hemisphere';
% allH(2).Title.String    = 'Right hemisphere';
% % allH(1).YTickLabel      = num2str((5:5:68)', '%02d');
% 
% % Overall legend
% ll = legend(allH(2), {'MATLAB', 'FEMA'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off');
% ll.Position(2) = ll.Position(2) - 0.11;
% ll.Position(1) = allH(2).Position(1) - 0.25;
% 
% % Save
% print('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/CorticalThickness_DK_FSE_nn/CorticalThickness_AgeDelta.png', '-dpng', '-r900');
% close(fig);
% 
% %% Make similar plot for age - cross sectional effect
% fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16], 'PaperUnits','centimeters', 'PaperSize', [16 16]);
% allH = tight_subplot(1, 2, [0.0 0.04], [0.12 0.03], [0.12 0.02]);
% 
% % Color scheme
% col_fitlme  = [217,95,2]./255;
% col_FEMA    = [27,158,119]./255;
% 
% % Locations to plot at
% xx1 = 0.7:1:34.7;
% xx2 = 1.3:1:35.3;
% 
% % Hold everything
% hold(allH(:), 'on');
% 
% % Settings
% markerType  = 'o';
% markerSize  = 3;
% LineStyle   = 'none';
% LineWidth   = 0.5;
% capSize     = 4;
% 
% alpha = 0.05/68;
% 
% % Plot
% % Slowly loop over all LH ROIs
% for rois = locLeft
%     % Determine if statistically significant or not - MATLAB
%     if pValues_MATLAB(locAge, locLeft(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(1), beta_MATLAB(locAge, locLeft(rois)), xx1(rois), beta_se_MATLAB(locAge, locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(1), beta_MATLAB(locAge, locLeft(rois)), xx1(rois), beta_se_MATLAB(locAge, locLeft(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% 
%     % Determine if statistically significant or not - FEMA
%     if pValues_FEMA(locAge, locLeft(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(1), beta_hat(locAge,    rois), xx2(rois), beta_se(locAge,        rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(1), beta_hat(locAge,    rois), xx2(rois), beta_se(locAge,        rois)*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% end
% 
% % Slowly loop over all RH rois
% for rois = 1:length(locRight)
%     % Determine if statistically significant or not - MATLAB
%     if pValues_MATLAB(locAge, locRight(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(2), beta_MATLAB(locAge, locRight(rois)), xx1(rois), beta_se_MATLAB(locAge, locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(2), beta_MATLAB(locAge, locRight(rois)), xx1(rois), beta_se_MATLAB(locAge, locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', 'k', 'MarkerFaceColor', 'none', 'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% 
%     % Determine if statistically significant or not - FEMA
%     if pValues_FEMA(locAge, locRight(rois)) < alpha
% 
%         % Plot with filled circle
%         errorbar(allH(2), beta_hat(locAge,    locRight(rois)), xx2(rois), beta_se(locAge,        locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
%     else
%         % Plot without filling
%         errorbar(allH(2), beta_hat(locAge,    locRight(rois)), xx2(rois), beta_se(locAge,        locRight(rois))*1.96, 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', 'k',   'MarkerFaceColor', 'none',   'LineWidth', LineWidth, 'CapSize', capSize);
%     end
% end
% 
% for ax = 1:2
%     % Apply axes settings
%     allH(ax).YLim            = [0 35];
%     allH(ax).Box             = 'off';
%     allH(ax).XLim            = [-0.026 0.026];
%     allH(ax).YTick           = 1:34;
%     allH(ax).YGrid           = 'on';
%     allH(ax).XLabel.String   = 'Estimated \beta age_{recruitment}';
%     allH(ax).Title.FontWeight = 'normal';
%     allH(ax).XAxis.FontSize  = 8;
%     allH(ax).XLabel.FontSize = 11;
%     allH(ax).YLabel.FontSize = 11;
%     allH(ax).Title.FontSize  = 11;
%     allH(ax).XTick           = -0.025:0.01:0.025;
%     allH(ax).XTickLabel      = -0.025:0.01:0.025;
% end
% 
% % allH(1).YLabel.String   = 'Imaging variables';
% allH(1).YTick           = 1:34;
% allH(1).YTickLabel      = strrep(roiNames_Use(locLeft), 'lh', '');
% % allH(1).TickLabelInterpreter = 'latex';
% 
% allH(1).Title.String    = 'Left hemisphere';
% allH(2).Title.String    = 'Right hemisphere';
% % allH(1).YTickLabel      = num2str((5:5:68)', '%02d');
% 
% % Overall legend
% ll = legend(allH(2), {'MATLAB', 'FEMA'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off');
% ll.Position(2) = ll.Position(2) - 0.11;
% ll.Position(1) = allH(2).Position(1) - 0.25;
% 
% % Save
% print('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/CorticalThickness_DK_FSE_nn/CorticalThickness_Age.png', '-dpng', '-r900');
% close(fig);
% 
% %% Plot random effects of family, subject, and error
% fig  = figure('Units', 'centimeters', 'Position', [10 10 16 22], 'PaperUnits','centimeters');
% allH = tight_subplot(1, 3, [0.0 0.01], [0.09 0.03], [0.11 0.01]);
% 
% % Color scheme
% col_fitlme  = [217,95,2]./255;
% col_FEMA    = [27,158,119]./255;
% 
% % Locations to plot at
% xx1 = 0.7:1:67.7;
% xx2 = 1.3:1:68.3;
% 
% % Hold all plots
% hold(allH(:), 'on');
% 
% % Compile lower and upper confidence interval for FEMA
% for params = 1:3
%     sig2Low(params,:) = prctile(squeeze(sig2mat_perm(params,:,2:end)) .* squeeze(sig2tvec_perm(:,:,2:end)),  2.5, 2);
%     sig2Upp(params,:) = prctile(squeeze(sig2mat_perm(params,:,2:end)) .* squeeze(sig2tvec_perm(:,:,2:end)), 97.5, 2);
% end
% 
% % Plot
% % Remember that MATLAB estimates are already squares
% for ax = 1:3
%     errorbar(allH(ax), varComp_MATLAB(ax,  1:end),      xx1, lowVar_MATLAB(ax, 1:end),  UppVar_MATLAB(ax, 1:end), 'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_fitlme, 'Color', col_fitlme, 'MarkerFaceColor', col_fitlme, 'LineWidth', LineWidth, 'CapSize', capSize);
%     errorbar(allH(ax), sig2mat(ax,   1:end).*sig2tvec,  xx2, sig2Low(ax,  1:end),       sig2Upp(ax,  1:end),      'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
%     % errorbar(allH(ax), sig2mat(ax,   1:end),       xx2, sig2Low(ax,  1:end),      sig2Upp(ax,  1:end),      'horizontal', 'LineStyle', LineStyle, 'Marker', markerType, 'MarkerSize', markerSize, 'MarkerEdgeColor', col_FEMA,   'Color', col_FEMA,   'MarkerFaceColor', col_FEMA,   'LineWidth', LineWidth, 'CapSize', capSize);
% 
%     % Apply axes settings
%     allH(ax).YLim               = [0 69];
%     allH(ax).Box                = 'off';
%     allH(ax).XLim               = [-0.1 1.7];
%     allH(ax).YTick              = 1:68;
%     allH(ax).YTickLabel         = [];
%     allH(ax).YAxis.TickLength   = [0 0];
%     allH(ax).YGrid              = 'on';
%     allH(ax).XTick              = 0:0.4:1.6;
%     allH(ax).XTickLabel         = 0:0.4:1.6;
%     allH(ax).YAxis.FontSize     = 8;
%     allH(ax).XAxis.FontSize     = 8;
%     allH(ax).YLabel.FontSize    = 11;
%     allH(ax).XLabel.FontSize    = 11;
%     allH(ax).Title.FontSize     = 11;
%     allH(ax).Title.FontWeight   = 'normal';
%     allH(ax).XLabel.String      = 'Estimated \sigma^2';
% end
% 
% allH(1).Title.String    = 'Family effect (\itF\rm)';
% allH(2).Title.String    = 'Subject effect (\itS\rm)';
% allH(3).Title.String    = 'Unmodeled variance (\itE\rm)';
% 
% allH(1).YTick           = 1:68;
% allH(1).YTickLabel       = strrep(strrep(roiNames_Use, 'lh', ' lh'), 'rh', ' rh');
% 
% % Overall legend
% ll              = legend(allH(2), {'MATLAB', 'FEMA'}, 'Orientation', 'horizontal', 'Location', 'southoutside', 'FontSize', 12, 'Box', 'off');
% ll.Position(2)  = ll.Position(2) - 0.09;
% 
% % Save
% print('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/CorticalThickness_DK_FSE_nn/CorticalThickness_RandomEffects.png', '-dpng', '-r900');
% close(fig);