%% Plot results for cortical thickness - vertex wise
addpath('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/cmig_tools-2.3.0/cmig_tools_utils/matlab/');
addpath('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/cmig_tools-2.3.0/showSurf');

% Settings
fname_results   = '/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/CorticalThickness_vertexWise_FSE/Results_masked.mat'; 

% Load results
load(fname_results);

% Load plotting template et al
load SurfView_surfs.mat % load surface templates
global icnum icnvert icsurfs curvvec_lh curvvec_rh surf_lh_pial surf_rh_pial %#ok<NUSED>
ico     = 5; % ico number
icnum   = ico+1; % index for ico number (icnum = ico + 1)
icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order

% Results to plot
locAge      = 2;
locAgeDelta = 3;

% Alpha value
numVertices = 18742;
alpha       = 0.05/numVertices;

% p values for age
pValues_age      = normcdf(-abs(zmat(locAge,:)))*2;
pValues_ageDelta = normcdf(-abs(zmat(locAgeDelta,:)))*2;

%% Age - uncorrected
% Define vertex values
vertvals = zmat(locAge, :);

fh = doPlot(vertvals, '\itZ\rm\bf-score age_{recruitment}', 'blueblackred', 2, -min(300, max(abs(vertvals))));
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/ZScores_Age_blueblackred_uncorrected.png', '-dpng', '-r900');
close(fh);

%% Age - Bonferroni corrected
% Define vertex values
vertvals = zmat(locAge, :);

% Set vertices which did not survive Bonferroni correction to zero
vertvals(~(pValues_age < alpha)) = 0;

fh = doPlot(vertvals, '\itZ\rm\bf-score age_{recruitment}', 'blueblackred', 2, -min(300, max(abs(vertvals))));
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/ZScores_Age_blueblackred_Bonferroni.png', '-dpng', '-r900');
close(fh);

%% AgeDelta - uncorrected
% Define vertex values
vertvals = zmat(locAgeDelta, :);
fh = doPlot(vertvals, '\itZ\rm\bf-score age_{delta}', 'blueblackred', 2, -min(300, max(abs(vertvals))));
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/ZScores_AgeDelta_blueblackred_uncorrected.png', '-dpng', '-r900');
close(fh);

%% AgeDelta - Bonferroni corrected
% Define vertex values
vertvals = zmat(locAgeDelta, :);

% Set vertices which did not survive Bonferroni correction to zero
vertvals(~(pValues_ageDelta < alpha)) = 0;

fh = doPlot(vertvals, '\itZ\rm\bf-score age_{delta}', 'blueblackred', 2, -min(300, max(abs(vertvals))));

print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/ZScores_AgeDelta_blueblackred_Bonferroni.png', '-dpng', '-r900');
close(fh);

%% Family effect = 1
fh = doPlot(sig2mat(1,:).*sig2tvec, 'Family Effect', 'fire', 1, 0);
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/RFX_FamilyEffect_fire.png', '-dpng', '-r900');
close(fh);

%% Subject effect = 2
fh = doPlot(sig2mat(2,:) .* sig2tvec, 'Subject Effect', 'fire', 1, 0);
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/RFX_SubjectEffect_fire.png', '-dpng', '-r900');
close(fh);

%% Error effect = 3
fh = doPlot(sig2mat(3,:) .* sig2tvec, 'Unmodeled variance', 'fire', 1, 0);
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/RFX_ErrorEffect_fire.png', '-dpng', '-r900');
close(fh);

%% Effect of F, and S
fh = doPlot(sum(sig2mat(1:2,:) .* sig2tvec), '\itF\rm\bf and \itS\rm\bf Effect', 'fire', 1, 0);
print('/Applications/Projects/2023-01-31_FEMAExperiments/2023-11-20_redone/Figures/CorticalThickness_vertexWise_FSE/RFX_F_S_Effect_fire.png', '-dpng', '-r900');
close(fh);

function fh = doPlot(vertvals, str, cmap, polarity, climMin)

% Inherit variables from global environment
global icnum icnvert icsurfs curvvec_lh curvvec_rh surf_lh_pial surf_rh_pial

vertvals_lh = vertvals(1:icnvert); % divide statistics by hemisphere for plotting
vertvals_rh = vertvals(icnvert+[1:icnvert]);

% specify limits for plot based on vertvals
fmax = min(300,max(abs(vertvals))); % max limit for plotting purposes
fmin = 0.0; % min limit for plotting purposes
fmid = fmax/2; % middle value for plotting purposes
fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function

curvcontrast = [0.2 0.2]; % contrast of gyri/sulci
bgcol        = [0 0 0]; % change to [1 1 1] for white background

% set colorbar limits - usually [fmin fmax] or [-fmax fmax]
clim = [climMin fmax];

% Define space
hvgap = [0.02 0.02];
lrgap = [0.02 0.02];
btgap = [0.12 0.01];

% Background - black
bgcol = [0 0 0];

% Create figure
fh = figure('Units', 'centimeters', 'Position', [10 10 16 10], 'Color', bgcol, 'InvertHardcopy', 'off');

% Create axes
allH = tight_subplot(2, 2, hvgap, btgap, lrgap);
hold(allH(:), 'on');

cm = eval(cmap);

axes(allH(1)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
axes(allH(2)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
axes(allH(3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
axes(allH(4)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;

% Set colorbar
colormap(cmap);
cb                  = colorbar('color', 'w');
cb.Label.String     = str; 
cb.Label.FontSize   = 12;
cb.FontSize         = 10;
cb.Box              = 'off';
cb.Label.FontWeight = 'bold';   
cb.Location         = 'south';
cb.Position(1)      = allH(1).Position(1);
cb.Position(2)      = cb.Position(2) - btgap(1);
cb.Position(3)      = allH(1).Position(3)*2 + hvgap(1);
caxis(clim);
end