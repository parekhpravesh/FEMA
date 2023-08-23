%% Plot corrmat results
load('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE/Results_half.mat', 'zmat');
load('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE/tmp_roinames.mat');
roinames    = roinames';
numParcels  = length(roinames);

%% Format results as a matrix
[r, c]  = find(tril(randn(numParcels, numParcels), -1));
idx     = sub2ind([numParcels numParcels], r, c);

zmat_age      = zeros(numParcels, numParcels);
zmat_ageDelta = zeros(numParcels, numParcels);

zmat_age(idx)       = zmat(2,:);
zmat_ageDelta(idx)  = zmat(3,:);

zmat_age        = triu(zmat_age.',1)        + tril(zmat_age);
zmat_ageDelta   = triu(zmat_ageDelta.',1)   + tril(zmat_ageDelta);

%% Gordon information
gordonNames     = readtable('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE/nodeNames.txt', 'ReadVariableNames', false);
gordonComm      = readtable('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE/communityAffiliation.txt', 'ReadVariableNames', false);
gordonCommNames = readtable('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE/communityNames.txt', 'ReadVariableNames', false);

% Parse Gordon parcellations
gordonNames.Var1      = strrep(strrep(gordonNames.Var1, 'L', 'Left '), 'R', 'Right ');
gordonNames.names     = strcat(gordonNames.Var1, gordonNames.Var2); 
gordonNames.Comm      = gordonComm.Var1;
gordonNames.CommNames = repmat({''}, height(gordonNames), 1);

for tmp = 1:height(gordonCommNames)
    loc = gordonNames.Comm == tmp;
    gordonNames.CommNames(loc) = gordonCommNames.Var1(tmp);
end

% Sort Gordon names by community affiliation
[gordonNames, rc] = sortrows(gordonNames, 'Comm');

% Some sort of order in ROI naming
% names = strrep(strrep(strrep(roinames, 'ctx-', ''), 'lh-', 'Left '), 'rh-', 'Right ');
% names = strrep(strrep(strrep(names,    'ctx_', ''), 'lh_', 'Left '), 'rh_', 'Right ');
% names = strrep(strrep(names,    'Left-', 'Left '), 'Right-', 'Right ');
% 
% loc_LH_DK = 1:34;
% loc_RH_DK = 35:68;
% loc_LR_DK = [loc_LH_DK; loc_RH_DK];
% loc_LR_DK = loc_LR_DK(:);
% 
% loc_LH_DT = 69:142;
% loc_RH_DT = 143:216;
% loc_LR_DT = [loc_LH_DT; loc_RH_DT];
% loc_LR_DT = loc_LR_DT(:);
% 
% loc_LH_GP = 217:377; % find(strcmpi(gordonNames.Var1, 'Left '));
% loc_RH_GP = 378:549; % find(strcmpi(gordonNames.Var1, 'Right '));
% loc_LR_GP = [loc_LH_GP, loc_RH_GP];
% loc_LR_GP = loc_LR_GP(rc);
% % loc_LR_GP = [loc_LH_GP; loc_RH_GP];
% % loc_LR_GP = loc_LR_GP(:);
% 
% loc_LH_aseg = [550:559, 563, 564, 566, 567];
% loc_RH_aseg = 568:581;
% loc_remain = [560:562, 565, 582];
% loc_LR_aseg = [loc_LH_aseg'; loc_RH_aseg'; loc_remain'];
% loc_LR_aseg = loc_LR_aseg(:);
% 
% loc_all = [loc_LR_DK; loc_LR_DT; loc_LR_GP'; loc_LR_aseg];

%% Atlas locations
loc_LR_DK   = 1:68;
loc_LR_DT   = 69:216;
loc_LR_GP   = 217:549;
loc_LR_GP   = loc_LR_GP(rc);
loc_LR_aseg = 550:582;
loc_all     = [loc_LR_DK, loc_LR_DT, loc_LR_GP, loc_LR_aseg];

% Gordon communities
loc_LR_GP_comm = cell(13, 1);
for comm = 1:13
    loc_LR_GP_comm{comm}  = loc_LR_GP((gordonNames.Comm == comm));
end

%% Number of parcels
np_DK       = length(loc_LR_DK);
np_DT       = length(loc_LR_DT);
np_GP       = length(loc_LR_GP);
np_aseg     = length(loc_LR_aseg);
np_GP_comm  = cellfun(@length, loc_LR_GP_comm);

%% Tick locations
% tickLoc = [1:2:np_DK, np_DK+1:2:np_DK+np_DT, np_DK+np_DT+1, np_DK+np_DT+np_GP+1:2:np_DK+np_DT+np_GP+np_aseg-length(loc_remain), np_DK+np_DT+np_GP+np_aseg-length(loc_remain)+1:numParcels];

%% Rearrange zmat
zmat_age        = zmat_age(loc_all,      loc_all);
zmat_ageDelta   = zmat_ageDelta(loc_all, loc_all);

%% Some settings
commWidth   = 1;
parcelWidth = 1.5;
parcelColor = [1 1 1];
commColor   = [0.75 0.75 0.75];

%% Locations where new parcellation starts
startPoint = [1,     np_DK+1,     np_DK+np_DT+1];
endPoint   = [np_DK, np_DK+np_DT];
for comm = 1:13
   startPoint(comm+3) = np_DK+np_DT+sum(np_GP_comm(1:comm))+1;
   endPoint(comm+2)   = np_DK+np_DT+sum(np_GP_comm(1:comm));
end
endPoint(end+1) = numParcels;
lenParcel  = endPoint - startPoint;

%% Age
limits = max(abs(zmat_age(:)));
fig    = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
imagesc(zmat_age);
colormap('blueblackred');
cm = colorbar('Location', 'eastoutside', 'Limits', [-limits, limits], 'Box', 'off');
cm.Label.String = '\itZ\rm-score age_{recruitment}';
cm.Label.FontSize = 10;
axis('equal');
set(gca, 'YDir', 'normal');
xlim([1 numParcels]);
ylim([1 numParcels]);

hold on

for lines = [1,2,16]
    rectangle('Position', [startPoint(lines), startPoint(lines), lenParcel(lines), lenParcel(lines)], 'LineWidth', parcelWidth, 'EdgeColor', parcelColor);
end
rectangle('Position', [np_DK+np_DT+1, np_DK+np_DT+1, np_GP-1, np_GP-1], 'LineWidth', parcelWidth, 'EdgeColor', parcelColor);
for lines = 3:15
    rectangle('Position', [startPoint(lines), startPoint(lines), lenParcel(lines), lenParcel(lines)], 'LineWidth', commWidth, 'EdgeColor', commColor);
end

xticks([]);
yticks([]);
ax = gca;
ax.XAxis.TickLength = [0 0];
ax.YAxis.TickLength = [0 0];
ax.XAxis.TickLabels = [];
ax.YAxis.TickLabels = [];

print('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE/Corrmat_Age_FSE.png', '-dpng', '-r900');
close(fig);

%% Age delta
limits = max(abs(zmat_ageDelta(:)));
fig    = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
imagesc(zmat_ageDelta);
colormap('blueblackred');
cm = colorbar('Location', 'eastoutside', 'Limits', [-limits, limits], 'Box', 'off');
cm.Label.String = '\itZ\rm-score age_{delta}';
cm.Label.FontSize = 10;
axis('equal');
set(gca, 'YDir', 'normal');
xlim([1 numParcels]);
ylim([1 numParcels]);

hold on

for lines = [1,2,16]
    rectangle('Position', [startPoint(lines), startPoint(lines), lenParcel(lines), lenParcel(lines)], 'LineWidth', parcelWidth, 'EdgeColor', parcelColor);
end
rectangle('Position', [np_DK+np_DT+1, np_DK+np_DT+1, np_GP-1, np_GP-1], 'LineWidth', parcelWidth, 'EdgeColor', parcelColor);
for lines = 3:15
    rectangle('Position', [startPoint(lines), startPoint(lines), lenParcel(lines), lenParcel(lines)], 'LineWidth', commWidth, 'EdgeColor', commColor);
end

xticks([]);
yticks([]);
ax = gca;
ax.XAxis.TickLength = [0 0];
ax.YAxis.TickLength = [0 0];
ax.XAxis.TickLabels = [];
ax.YAxis.TickLabels = [];

print('/home/pparekh/analyses/2023-02-17_FEMA-ABCD/Corrmat_FSE/Corrmat_AgeDelta_FSE.png', '-dpng', '-r900');
close(fig);

% %% Age
% limits = max(abs(zmat_age(:)));
% fig    = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
% imagesc(zmat_age(loc_all, loc_all));
% colormap('blueblackred');
% colorbar('Location', 'eastoutside', 'Limits', [-limits, limits]);
% axis('equal');
% set(gca, 'YDir', 'normal');
% xlim([1 numParcels]);
% ylim([1 numParcels]);
% 
% hold on
% % Draw boxes around parcellation schemes
% % DK
% plot(min(loc_LR_DK):max(loc_LR_DK),     repmat(min(loc_LR_DK), np_DK, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(min(loc_LR_DK), np_DK, 1),  min(loc_LR_DK):max(loc_LR_DK),      'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(min(loc_LR_DK):max(loc_LR_DK),     repmat(max(loc_LR_DK), np_DK, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(max(loc_LR_DK), np_DK, 1),  min(loc_LR_DK):max(loc_LR_DK),      'LineWidth', parcelWidth, 'Color', parcelColor);
% 
% % DT
% plot(min(loc_LR_DT):max(loc_LR_DT),     repmat(min(loc_LR_DT), np_DT, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(min(loc_LR_DT), np_DT, 1),  min(loc_LR_DT):max(loc_LR_DT),      'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(min(loc_LR_DT):max(loc_LR_DT),     repmat(max(loc_LR_DT), np_DT, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(max(loc_LR_DT), np_DT, 1),  min(loc_LR_DT):max(loc_LR_DT),      'LineWidth', parcelWidth, 'Color', parcelColor);
% 
% % GP
% plot(min(loc_LR_GP):max(loc_LR_GP),     repmat(min(loc_LR_GP), np_GP, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(min(loc_LR_GP), np_GP, 1),  min(loc_LR_GP):max(loc_LR_GP),      'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(min(loc_LR_GP):max(loc_LR_GP),     repmat(max(loc_LR_GP), np_GP, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(max(loc_LR_GP), np_GP, 1),  min(loc_LR_GP):max(loc_LR_GP),      'LineWidth', parcelWidth, 'Color', parcelColor);
% 
% % aseg
% plot(min(loc_LR_aseg):max(loc_LR_aseg),     repmat(min(loc_LR_aseg), np_aseg, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(min(loc_LR_aseg), np_aseg, 1),  min(loc_LR_aseg):max(loc_LR_aseg),      'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(min(loc_LR_aseg):max(loc_LR_aseg),     repmat(max(loc_LR_aseg), np_aseg, 1),   'LineWidth', parcelWidth, 'Color', parcelColor);
% plot(repmat(max(loc_LR_aseg), np_aseg, 1),  min(loc_LR_aseg):max(loc_LR_aseg),      'LineWidth', parcelWidth, 'Color', parcelColor);
% 
% % Lines for Gordon parcels
% for comm = 1:13
%     plot([min(loc_LR_GP_comm{comm}) min(loc_LR_GP_comm{comm})+np_GP_comm(comm)],     [min(loc_LR_GP_comm{comm}), min(loc_LR_GP_comm{comm})],   'LineWidth', parcelWidth, 'Color', parcelColor);
%     plot([min(loc_LR_GP_comm{comm}) min(loc_LR_GP_comm{comm})],     [min(loc_LR_GP_comm{comm})  min(loc_LR_GP_comm{comm})+np_GP_comm(comm)],'LineWidth', parcelWidth, 'Color', parcelColor);
%     plot([min(loc_LR_GP_comm{comm})+np_GP_comm(comm) min(loc_LR_GP_comm{comm})+np_GP_comm(comm)],     [min(loc_LR_GP_comm{comm})  min(loc_LR_GP_comm{comm})+np_GP_comm(comm)],   'LineWidth', parcelWidth, 'Color', parcelColor);
%     plot([min(loc_LR_GP_comm{comm}) min(loc_LR_GP_comm{comm})+np_GP_comm(comm)],  [min(loc_LR_GP_comm{comm})+np_GP_comm(comm) min(loc_LR_GP_comm{comm})+np_GP_comm(comm)],      'LineWidth', parcelWidth, 'Color', parcelColor);
% end
% 
% % rectangle('Position', [min(loc_LR_DK)-0.5,   min(loc_LR_DK)-0.5,    max(loc_LR_DK)-min(loc_LR_DK)+1,      max(loc_LR_DK)-min(loc_LR_DK)+1],     'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% % rectangle('Position', [min(loc_LR_DT)-0.5,   min(loc_LR_DT)-0.5,    max(loc_LR_DT)-min(loc_LR_DT)+1,      max(loc_LR_DT)-min(loc_LR_DT)+1],     'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% % rectangle('Position', [min(loc_LR_GP)-0.5,   min(loc_LR_GP)-0.5,    max(loc_LR_GP)-min(loc_LR_GP)+1,      max(loc_LR_GP)-min(loc_LR_GP)+1],     'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% % rectangle('Position', [min(loc_LR_aseg)-0.5, min(loc_LR_aseg)-0.5,  max(loc_LR_aseg)-min(loc_LR_aseg)+1,  max(loc_LR_aseg)-min(loc_LR_aseg)+1], 'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% 
% % % Draw Gordon lines
% % count = 217;
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_1)+0.5,       length(loc_LR_GP_1)]+0.5,      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_1);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_2)+0.5,       length(loc_LR_GP_2)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_2);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_3)+0.5,       length(loc_LR_GP_3)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_3);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_4)+0.5,       length(loc_LR_GP_4)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_4);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_5)+0.5,       length(loc_LR_GP_5)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_5);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_6)+0.5,       length(loc_LR_GP_6)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_6);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_7)+0.5,       length(loc_LR_GP_7)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_7);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_8)+0.5,       length(loc_LR_GP_8)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_8);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_9)+0.5,       length(loc_LR_GP_9)+0.5],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_9);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_10)+0.5,      length(loc_LR_GP_10)+0.5],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_10);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_11)+0.5,      length(loc_LR_GP_11)+0.5],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_11);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_12)+0.5,      length(loc_LR_GP_12)+0.5],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % count = count + length(loc_LR_GP_12);
% % rectangle('Position', [count-0.5,    count-0.5,     length(loc_LR_GP_13)+0.5,      length(loc_LR_GP_13)+0.5],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% % 
% % xticks(sort([min(loc_LR_DK); min(loc_LR_DT); min(loc_LR_GP'); min(loc_LR_aseg); max(loc_LR_DK); max(loc_LR_DT); max(loc_LR_GP'); max(loc_LR_aseg)]));
% 
% %% Age delta
% limits = max(abs(zmat_ageDelta(:)));
% fig    = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
% imagesc(zmat_ageDelta(loc_all, loc_all));
% colormap('blueblackred');
% colorbar('Location', 'eastoutside', 'Limits', [-limits, limits]);
% axis('equal');
% set(gca, 'YDir', 'normal');
% xlim([1 numParcels]);
% ylim([1 numParcels]);
% 
% % Draw boxes around parcellation schemes
% rectangle('Position', [min(loc_LR_DK),   min(loc_LR_DK),    max(loc_LR_DK)-min(loc_LR_DK)+1,      max(loc_LR_DK)-min(loc_LR_DK)+1],     'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% rectangle('Position', [min(loc_LR_DT),   min(loc_LR_DT),    max(loc_LR_DT)-min(loc_LR_DT)+1,      max(loc_LR_DT)-min(loc_LR_DT)+1],     'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% rectangle('Position', [min(loc_LR_GP),   min(loc_LR_GP),    max(loc_LR_GP)-min(loc_LR_GP)+1,      max(loc_LR_GP)-min(loc_LR_GP)+1],     'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% rectangle('Position', [min(loc_LR_aseg), min(loc_LR_aseg),  max(loc_LR_aseg)-min(loc_LR_aseg)+1,  max(loc_LR_aseg)-min(loc_LR_aseg)+1], 'EdgeColor', [1, 1 ,1], 'LineWidth', 2);
% 
% % Draw Gordon lines
% count = 217;
% rectangle('Position', [count,    count,     length(loc_LR_GP_1),       length(loc_LR_GP_1)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_1);
% rectangle('Position', [count,    count,     length(loc_LR_GP_2),       length(loc_LR_GP_2)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_2);
% rectangle('Position', [count,    count,     length(loc_LR_GP_3),       length(loc_LR_GP_3)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_3);
% rectangle('Position', [count,    count,     length(loc_LR_GP_4),       length(loc_LR_GP_4)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_4);
% rectangle('Position', [count,    count,     length(loc_LR_GP_5),       length(loc_LR_GP_5)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_5);
% rectangle('Position', [count,    count,     length(loc_LR_GP_6),       length(loc_LR_GP_6)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_6);
% rectangle('Position', [count,    count,     length(loc_LR_GP_7),       length(loc_LR_GP_7)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_7);
% rectangle('Position', [count,    count,     length(loc_LR_GP_8),       length(loc_LR_GP_8)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_8);
% rectangle('Position', [count,    count,     length(loc_LR_GP_9),       length(loc_LR_GP_9)],      'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_9);
% rectangle('Position', [count,    count,     length(loc_LR_GP_10),      length(loc_LR_GP_10)],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_10);
% rectangle('Position', [count,    count,     length(loc_LR_GP_11),      length(loc_LR_GP_11)],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_11);
% rectangle('Position', [count,    count,     length(loc_LR_GP_12),      length(loc_LR_GP_12)],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% count = count + length(loc_LR_GP_12);
% rectangle('Position', [count,    count,     length(loc_LR_GP_13),      length(loc_LR_GP_13)],     'EdgeColor', [.5, .5 ,.5], 'LineWidth', 1);
% %%
% fig  = figure('Units', 'centimeters', 'Position', [10 10 16 16]);
% allH = tight_subplot(1, 1, 0.05, [0.02 0.01], [0.08 0.08]);
% imagesc(allH(1), zmat_age); % .* (pValues_age(271:549, 271:549) < 0.05));
% colormap(allH(1), 'blueblackred');
% colorbar(allH(1), 'Location', 'eastoutside', 'Limits', [-max(abs(zmat_age(:))), max(abs(zmat_age(:)))]);
% axis(allH(1), 'equal');
% allH(1).YDir = 'normal';
% xlim(allH(1), [1 numParcels]);
% ylim(allH(1), [1 numParcels]);
