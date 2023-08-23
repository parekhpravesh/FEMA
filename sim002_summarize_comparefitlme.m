%% Summarize results of simulation 002: comparison with fitlme

% Set paths and initialize
rootDir         = '/ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/Simulation002_Compare_fitlme_perms_nn';
settings        = load(fullfile(rootDir, 'Settings.mat'), 'settings');
settings        = settings.settings;
results_fitlme  = cell(10, 1);
results_FEMA    = cell(10, 1);

% Go over every number of observation and extract parameters
for vals            = 1:10
    toLoad          = fullfile(rootDir, ['nObs-', num2str(settings.nObservations(vals), '%02d')]);
    vars_FEMA       = load(fullfile(toLoad, 'Results_FEMA.mat'), 'beta_*', 'FEMAelapsed', 'sig2*');
    vars_fitlme     = load(fullfile(toLoad, 'Results_fitlme.mat'), 'mdls', 'lmeElapsed');
    
    % Extract out random effects
    [beta_hat_fitlme, beta_se_fitlme]                   = deal(zeros(settings.nXvars, settings.nyVars));
    [sig2mat_fitlme,  sig2Low_fitlme, sig2Upp_fitlme]   = deal(zeros(3, settings.nyVars));
    
    for mdls = 1:settings.nyVars
        [psi, mse, stats]                        = covarianceParameters(vars_fitlme.mdls{mdls,1});
        beta_hat_fitlme(1:settings.nXvars, mdls) = vars_fitlme.mdls{mdls,1}.Coefficients.Estimate;
        beta_se_fitlme(1:settings.nXvars,  mdls) = vars_fitlme.mdls{mdls,1}.Coefficients.SE;
        sig2mat_fitlme(1:3,                mdls) = [stats{1}.Estimate; stats{2}.Estimate; stats{3}.Estimate];
        sig2Low_fitlme(1:3,                mdls) = [stats{1}.Lower; stats{2}.Lower; stats{3}.Lower];
        sig2Upp_fitlme(1:3,                mdls) = [stats{1}.Upper; stats{2}.Upper; stats{3}.Upper];
    end
    
    % Remember that these are in standard deviation units - to be squared
    % before comparing with FEMA output
    variables_lme.beta_hat_fitlme = beta_hat_fitlme;
    variables_lme.beta_se_fitlme  = beta_se_fitlme;
    variables_lme.sig2mat_fitlme  = sig2mat_fitlme;
    variables_lme.sig2Low_fitlme  = sig2Low_fitlme;
    variables_lme.sig2Upp_fitlme  = sig2Upp_fitlme;
    variables_lme.elapsed         = vars_fitlme.lmeElapsed;
    results_fitlme{vals, 1}       = variables_lme;
    
    % Compile for FEMA - remember to convert back to real scale by
    % multiplying with sig2tvec
    variables_FEMA.beta_hat_FEMA    = vars_FEMA.beta_hat(2:end, :);
    variables_FEMA.beta_se_FEMA     = vars_FEMA.beta_se(2:end,  :);
    variables_FEMA.sig2mat_FEMA     = vars_FEMA.sig2mat .* vars_FEMA.sig2tvec;
    
    for params  = 1:3
        variables_FEMA.sig2Low_FEMA(params,:) = prctile(squeeze(vars_FEMA.sig2mat_perm(params,:,2:end)) .* squeeze(vars_FEMA.sig2tvec_perm(:,:,2:end)),  2.5, 2);
        variables_FEMA.sig2Upp_FEMA(params,:) = prctile(squeeze(vars_FEMA.sig2mat_perm(params,:,2:end)) .* squeeze(vars_FEMA.sig2tvec_perm(:,:,2:end)), 97.5, 2);
    end
    variables_FEMA.elapsed  = vars_FEMA.FEMAelapsed;
    results_FEMA{vals, 1}   = variables_FEMA;
end

save(fullfile(rootDir, 'Results_experiment02_nn.mat'), 'results_FEMA', 'results_fitlme', 'settings');