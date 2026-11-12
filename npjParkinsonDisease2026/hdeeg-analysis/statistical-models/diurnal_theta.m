% =========================================================================
% EEG ANALYSIS: LINEAR MIXED-EFFECTS MODEL WITH GROUP × SESSION INTERACTION
%
% This script performs statistical analysis on EEG theta-band power data
% using linear mixed-effects models (LME) to compare a reference group 
% (e.g., CTL) against clinical groups (EPD, ADV, DYS).
%
% Analysis includes:
% - Group comparisons (CTL vs ADV, EPD, DYS)
% - Group × Session interaction (morning vs evening)
% - Integration of clinical covariates (age, sex)
% - Calculation of effect sizes (Cohen’s d) based on total variance
%
% Variables:
% - Theta-band PSD computed for each EEG channel
% - Clinical covariates loaded from an Excel file
% - EEG channels filtered to exclude noisy electrodes
%
% Output:
% - t-values, p-values, and Cohen’s d for:
%     - Group × Session interaction
%     - Main Group effect
%     - Main Session effect
% - mEEG (marginal EEG) and design matrix X for downstream analysis
% - Result structure saved for each group comparison
%
% Requirements:
% - EEGLAB (for .ced channel location file)
% - PSD file structured by subject/group/session
% - Excel file with subject IDs and clinical variables
%
% Author: [Giovanni Lombardi]
% Institution: [Institute of Clinical Neuroscience of Southern Switzerland]
% Date: December[2025]
%
% Statistical significance of the Group × Session interaction effect 
% was assessed using a non-parametric permutation test with 
% Threshold-Free Cluster Enhancement (TFCE).
% For each EEG channel, t-values from LME models were transformed via TFCE 
% and compared against a null distribution generated from 2000 permutations 
% of the group and session labels. p-values were computed as the proportion 
% of permutations where TFCE values exceeded the observed value, with 
% multiple comparisons correction based on the distribution of the 
% maximum TFCE values (Smith & Nichols, 2009). Channels above the 95th 
% percentile threshold were considered statistically significant.
% =========================================================================


%% [1] Load PSD data
load('...\psdIndividual.mat');  % psdIndividual structure

%% [2] Setup: Reference group and session
referenceGroup = 'CTL';         % Reference group
referenceSession = 'morning';   % Reference session
fprintf('Selected baseline: Group = %s | Session = %s\n', referenceGroup, referenceSession);

% Load and filter channel locations
chanlocs_all = readlocs('...\ANT_64.ced');
exclude = [13 19 32];  % Remove noisy or irrelevant channels
chanlocs = chanlocs_all;
chanlocs(exclude) = [];

%% [3] Build PSD array from structure
groupNames = {'CTL', 'EPD', 'ADV', 'DYS'};
[psdArray, labels_study] = build_PSD_data_from_struct(psdIndividual, groupNames);
psdArray = log10(psdArray);  % Convert to dB scale

% Extract study metadata
Id      = labels_study.Id;
Group   = labels_study.Group;
Session = labels_study.Session;

% Extract theta-band power (4–8 Hz)
theta_idx = find(freqs >= 4 & freqs <= 8);
fprintf('Theta frequency range: %.2f – %.2f Hz\n', freqs(theta_idx(1)), freqs(theta_idx(end)));
psd_theta = squeeze(mean(psdArray(:, theta_idx, :), 2));  % Resulting size: [channels × subjects]
nCh = size(psd_theta, 1);

%% === STEP 4 – Fit Group × Session model
compareList = {'EPD','ADV','DYS'};
%compareList = {'CTL','EPD','ADV'};

results_interaction = struct;

% === Load clinical data
fileClin = '...\clinical_data_old.xlsx';
Tclin = readtable(fileClin, 'VariableNamingRule', 'preserve');
Tclin.ID = replace(lower(string(Tclin.ID)), 'ta', 'x');  % normalize IDs

for g = 1:length(compareList)
    groupLabel = compareList{g};
    fprintf('\n>>> Comparison %s vs %s\n', referenceGroup, groupLabel);

    idx = Group == referenceGroup | Group == groupLabel;
    Id_cmp      = Id(idx);
    Group_cmp   = Group(idx);
    Session_cmp = Session(idx);
    psd_cmp     = psd_theta(:, idx);

    tvals = nan(nCh,1); cohen_d = nan(nCh,1);
    tvals_group = nan(nCh,1); cohen_d_group = nan(nCh,1);
    tvals_session = nan(nCh,1); cohen_d_session = nan(nCh,1);
    pval_interaction = nan(nCh,1); pval_group = nan(nCh,1); pval_session = nan(nCh,1);
    betas_all = nan(nCh, 6); se_all = nan(nCh, 6);
    mEEG_all = nan(nCh, length(Id_cmp)); X_ref = [];

    for ch = 1:nCh
        y = double(psd_cmp(ch, :)');
        T = table(y, Id_cmp, Session_cmp, Group_cmp, 'VariableNames', {'y','Id','Session','Group'});
        T.Id = string(T.Id);
        T.Session = categorical(string(T.Session), {'morning','evening'});
        T.Session = reordercats(T.Session, [referenceSession, setdiff({'morning','evening'}, referenceSession)]);
        T.Group = categorical(string(T.Group), {referenceGroup, groupLabel});

        id_clean = lower(replace(string(T.Id), 'ta', 'x'));
        [found, idxC] = ismember(id_clean, Tclin.ID);

        T.age = NaN(height(T),1);
        T.sex = NaN(height(T),1);
        T.age(found) = Tclin.Age(idxC(found));
        T.sex(found) = double(strcmpi(Tclin.Sex(idxC(found)), 'M'));

        try
            mdl = fitlme(T, 'y ~ Group*Session + age + sex + (1|Id)');
            X = designMatrix(mdl);
            beta_row = fixedEffects(mdl);
            t_row = mdl.Coefficients.tStat;
            p_values = mdl.Coefficients.pValue;
            mEEG = fitted(mdl, 'Conditional', false) + residuals(mdl);

            if isempty(X_ref)
                X_ref = X;
            elseif size(X_ref,2) ~= size(X,2)
                warning("Channel %d – design matrix X has a different size than X_ref", ch);
            end

            otherSession = setdiff({'morning','evening'}, referenceSession);
           
            if strcmp(referenceSession, 'morning')
                    coeffNameA = sprintf('Group_%s:Session_evening', groupLabel);
                    coeffNameB = sprintf('Session_evening:Group_%s', groupLabel);
            else
                    coeffNameA = sprintf('Group_%s:Session_morning', groupLabel);
                    coeffNameB = sprintf('Session_morning:Group_%s', groupLabel);
            end

            idx_inter = find(strcmp(mdl.CoefficientNames, coeffNameA) | strcmp(mdl.CoefficientNames, coeffNameB));

            var_intercept_part = 0; var_slope_part = 0;
            C = mdl.covarianceParameters;
            if iscell(C)
                for k = 1:length(C)
                    if isstruct(C{k})
                        var_intercept_part = var_intercept_part + C{k}(1).Estimate^2;
                        if numel(C{k}) > 1
                            var_slope_part = var_slope_part + C{k}(2).Estimate^2;
                        end
                    elseif isfield(C{k}, 'Estimate')
                        var_intercept_part = var_intercept_part + C{k}.Estimate^2;
                    end
                end
            elseif isstruct(C)
                var_intercept_part = C(1).Estimate^2;
                if numel(C) > 1
                    var_slope_part = C(2).Estimate^2;
                end
            end

            var_residual = mdl.MSE;
            denom = sqrt(var_intercept_part + var_slope_part + var_residual);

            if ~isempty(idx_inter)
                tvals(ch) = t_row(idx_inter);
                pval_interaction(ch) = p_values(idx_inter);
                cohen_d(ch) = beta_row(idx_inter) / denom;
            end

            coeffGroup = sprintf('Group_%s', groupLabel);
            idx_group = find(strcmp(mdl.CoefficientNames, coeffGroup));
            if ~isempty(idx_group)
                tvals_group(ch) = t_row(idx_group);
                pval_group(ch) = p_values(idx_group);
                cohen_d_group(ch) = beta_row(idx_group) / denom;
            end

            % === Main effect of Session
            if strcmp(referenceSession, 'morning')
                sessionCoeff = 'Session_evening';
            else
                sessionCoeff = 'Session_morning';
            end

            idx_session = find(strcmp(mdl.CoefficientNames, sessionCoeff));
            if ~isempty(idx_session)
                tvals_session(ch) = t_row(idx_session);
                pval_session(ch) = p_values(idx_session);
                cohen_d_session(ch) = beta_row(idx_session) / denom;
            else
                tvals_session(ch) = NaN;
                pval_session(ch) = NaN;
                cohen_d_session(ch) = 0;
            end

            betas_all(ch, 1:length(beta_row)) = beta_row;
            se_all(ch, 1:length(t_row)) = mdl.Coefficients.SE;
            mEEG_all(ch,:) = mEEG';

        catch ME
            warning("Error in %s vs reference group – Channel %d: %s", groupLabel, ch, ME.message);
        end
    end

    results_interaction.(groupLabel).tvals = tvals;
    results_interaction.(groupLabel).cohen_d = cohen_d;
    results_interaction.(groupLabel).tvals_group = tvals_group;
    results_interaction.(groupLabel).cohen_d_group = cohen_d_group;
    results_interaction.(groupLabel).tvals_session = tvals_session;
    results_interaction.(groupLabel).cohen_d_session = cohen_d_session;
    results_interaction.(groupLabel).pval_session = pval_session;
    results_interaction.(groupLabel).pval_interaction = pval_interaction;
    results_interaction.(groupLabel).pval_group = pval_group;
    results_interaction.(groupLabel).betas = betas_all;
    results_interaction.(groupLabel).se = se_all;
    results_interaction.(groupLabel).mEEG = mEEG_all;
    results_interaction.(groupLabel).X = X_ref;
    results_interaction.(groupLabel).Id = Id_cmp;
    results_interaction.(groupLabel).coeffName = sprintf('Interaction: %s vs %s', referenceGroup, groupLabel);
end


%% === SETUP FOR PERMUTATIONS + TFCE === group*session
nPerms = 2000;
E_H = [0.66 2];
tfce_results_interaction = struct;
referenceGroup = 'CTL';

progressStep = round(nPerms / 20);
fprintf('TFCE Permutations (%d)...\n', nPerms);

for g = 1:length(compareList)
    groupLabel = compareList{g};
    tag = sprintf('TFCE_%s_vs_%s', referenceGroup, groupLabel);
    fprintf('\n%s\n', tag);

    idx = Group == referenceGroup | Group == groupLabel;
    Id_cmp      = Id(idx);
    Group_cmp   = Group(idx);
    Session_cmp = Session(idx);
    psd_cmp     = psd_theta(:, idx);
    nSubj       = sum(idx);
    nCh         = size(psd_cmp, 1);

   % Identify unique subjects
    [uniqueIDs, ~, subjectIdx] = unique(string(Id_cmp), 'stable');
    nSubjects = numel(uniqueIDs);

    % Extract one group label per subject
    subjectGroups = strings(nSubjects, 1);

    for s = 1:nSubjects
        gLabels = unique(string(Group_cmp(subjectIdx == s)));
        assert(numel(gLabels) == 1, 'Inconsistent group labels for subject');
        subjectGroups(s) = gLabels;
    end

    % === Extract clinical covariates
    id_clean = lower(replace(string(Id_cmp), 'ta', 'x'));
    [found, idxC] = ismember(id_clean, Tclin.ID);
    age_all = NaN(nSubj,1);
    sex_all = NaN(nSubj,1);
    age_all(found) = Tclin.Age(idxC(found));
    sex_all(found) = double(strcmpi(Tclin.Sex(idxC(found)), 'M'));

    % === Compute T_Obs (observed t-statistics)
    t_obs = nan(nCh, 1);
    for ch = 1:nCh
        y = double(psd_cmp(ch, :)');
        T = table(y, Id_cmp, Session_cmp, Group_cmp, age_all, sex_all, ...
            'VariableNames', {'y','Id','Session','Group','age','sex'});
        T.Id = string(T.Id);
        T.Session = categorical(string(T.Session), {'morning','evening'});
        T.Session = reordercats(T.Session, [referenceSession, setdiff({'morning','evening'}, referenceSession)]);
        T.Group = categorical(string(T.Group), {referenceGroup, groupLabel});

        if strcmp(referenceSession, 'morning')
            coeffNameA = sprintf('Group_%s:Session_evening', groupLabel);
            coeffNameB = sprintf('Session_evening:Group_%s', groupLabel);
        else
            coeffNameA = sprintf('Group_%s:Session_morning', groupLabel);
            coeffNameB = sprintf('Session_morning:Group_%s', groupLabel);
        end

        try
            mdl = fitlme(T, 'y ~ Group*Session + age + sex + (1|Id)');
            idx_inter = find(strcmp(mdl.CoefficientNames, coeffNameA) | strcmp(mdl.CoefficientNames, coeffNameB));
            if ~isempty(idx_inter)
                t_obs(ch) = mdl.Coefficients.tStat(idx_inter);
            end
        catch
        end
    end

    % === Permutations
    t_perms = nan(nPerms, nCh);
    parfor p = 1:nPerms
        t_perm_ch = nan(nCh, 1);
        shuff_idx = randperm(nSubjects);
        shuffledGroups = subjectGroups(shuff_idx);
        Group_perm = shuffledGroups(subjectIdx);
        Session_perm = Session_cmp;

        for ch = 1:nCh
            y = double(psd_cmp(ch, :)');
            T = table(y, Id_cmp, Session_perm, Group_perm, age_all, sex_all, ...
                'VariableNames', {'y','Id','Session','Group','age','sex'});
            T.Id = string(T.Id);
            T.Session = categorical(string(T.Session), {'morning','evening'});
            T.Session = reordercats(T.Session, [referenceSession, setdiff({'morning','evening'}, referenceSession)]);
            T.Group = categorical(string(T.Group), {referenceGroup, groupLabel});

            if strcmp(referenceSession, 'morning')
                coeffNameA = sprintf('Group_%s:Session_evening', groupLabel);
                coeffNameB = sprintf('Session_evening:Group_%s', groupLabel);
            else
                coeffNameA = sprintf('Group_%s:Session_morning', groupLabel);
                coeffNameB = sprintf('Session_morning:Group_%s', groupLabel);
            end

            try
                mdl = fitlme(T, 'y ~ Group*Session + age + sex + (1|Id)');
                idx_inter = find(strcmp(mdl.CoefficientNames, coeffNameA) | strcmp(mdl.CoefficientNames, coeffNameB));
                if ~isempty(idx_inter)
                    t_perm_ch(ch) = mdl.Coefficients.tStat(idx_inter);
                end
            catch
            end
        end

        t_perms(p,:) = t_perm_ch';
        if mod(p, progressStep) == 0
            fprintf('%d/%d permutations completed...\n', p, nPerms);
        end
    end

    Results = lmeEEG_TFCE_mod(t_obs, t_perms, chanlocs, E_H);
    cohen_d = results_interaction.(groupLabel).cohen_d;

    tfce_results_interaction.(tag).T_Obs      = t_obs;
    tfce_results_interaction.(tag).TFCE_Obs   = Results.TFCE_Obs;
    tfce_results_interaction.(tag).P_Values   = Results.P_Values;
    tfce_results_interaction.(tag).Mask       = Results.Mask;
    tfce_results_interaction.(tag).Cohen_d    = cohen_d;
    tfce_results_interaction.(tag).TFCE_Perms = Results.TFCE_Perms;
    tfce_results_interaction.(tag).maxTFCE    = Results.maxTFCE;
end


%% === PLOT TFCE PERMUTATION RESULTS ===

% Set the order of group comparisons (by disease severity)
compareList = {'EPD', 'ADV', 'DYS'};

% === Determine common color scale
max_TFCE = -inf;
for g = 1:length(compareList)
    tag = sprintf('TFCE_%s_vs_%s', referenceGroup, compareList{g});
    res = tfce_results_interaction.(tag);
    max_TFCE = max(max_TFCE, max(abs(res.TFCE_Obs(:))));
end

% Display labels of significant channels
disp(strjoin({chanlocs(find(res.Mask)).labels}, ', '));

%% === Generate 2x3 topoplot layout with titles ===
figure('Name', 'TFCE Topoplots', 'Color','w', 'Position',[100 100 1200 700]);
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(sprintf('TFCE – Group × Session Interaction\nBaseline: %s – %s', referenceGroup, referenceSession), 'FontSize', 14, 'FontWeight','bold');

for g = 1:length(compareList)
    groupLabel = compareList{g};
    tag = sprintf('TFCE_%s_vs_%s', referenceGroup, groupLabel);
    res = tfce_results_interaction.(tag);

    % Row 1: TFCE statistics
    nexttile(g);
    topoplot(res.TFCE_Obs, chanlocs, ...
        'maplimits', [-max_TFCE, max_TFCE], ...
        'electrodes', 'on', ...
        'style', 'both');
    colorbar;
    title(sprintf('TFCE: %s vs %s', groupLabel, referenceGroup), 'FontSize', 10);

    % Row 2: -log10(p) with significance markers
    nexttile(g + 3);
    pvals_log = -log10(res.P_Values);
    sigCh = find(res.Mask);
    topoplot(pvals_log, chanlocs, ...
        'maplimits', [0, 3], ...
        'electrodes', 'on', ...
        'style', 'both', ...
        'emarker2', {sigCh, 'o', [1 0 0], 4, 2});
    colorbar;
    title('-log10(p)', 'FontSize', 10);
end

% Export high-resolution figure
exportgraphics(gcf, 'Brain_clean_Publication_Figure.png', 'Resolution', 600);
