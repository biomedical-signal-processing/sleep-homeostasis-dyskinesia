% =========================================================================
% EEG ANALYSIS: BETWEEN-GROUP COMPARISON (MORNING SESSION ONLY)
%
% This script performs a statistical analysis on EEG theta power data,
% comparing the control group (CTL) with each clinical group (EPD, ADV, DYS)
% during the **morning session only**, using:
%
% 1. A linear model (fitlm) including the following covariates:
%    - Age
%    - Sex
%
% 2. A permutation-based inference procedure (n = 2000) for statistical validation,
%    followed by correction using TFCE (Threshold-Free Cluster Enhancement).
%
% Outputs:
% - Observed t-values for each channel
% - p-values
% - Cohen’s d effect sizes
% - TFCE statistics and significant channel masks
% - 2×3 topoplot visualization (TFCE and -log10(p))
%
% Required files:
% - `psdIndividual.mat` containing theta power data per subject
% - `clinical_data.xlsx` with subject ID, age, and sex
% - `ANT_64.ced` for electrode position coordinates
% Author: [Giovanni Lombardi]
% Institution: [Institute of Clinical Neuroscience of Southern Switzerland]
% Date: December[2025]
% =========================================================================

load('...\psdIndividual.mat');  % structure psdIndividual

%% === INITIAL SETUP
referenceGroup = 'CTL';        % Reference group
referenceSession = 'morning';  % Reference session

fprintf('Selected baseline: Group = %s | Session = %s\n', referenceGroup, referenceSession);

%% === STEP 0: Load and filter channels
chanlocs_all = readlocs('...\ANT_64.ced');
exclude = [13 19 32];
chanlocs = chanlocs_all;
chanlocs(exclude) = [];

%% === STEP 1: Load data
groupNames = {'CTL', 'EPD', 'ADV', 'DYS'};
[psdArray, labels_study] = build_PSD_data_from_struct(psdIndividual, groupNames);
psdArray = log10(psdArray);  % Convert to dB

%% === STEP 2: Extract study information
Id      = labels_study.Id;
Group   = labels_study.Group;
Session = labels_study.Session;

%% === STEP 3: Compute theta power (4–8 Hz)
theta_idx = find(freqs >= 4 & freqs <= 8);
fprintf('Theta frequencies used: %.2f - %.2f Hz\n', freqs(theta_idx(1)), freqs(theta_idx(end)));

psd_theta = squeeze(mean(psdArray(:, theta_idx, :), 2));  % [channels x subjects]
nCh = size(psd_theta, 1);

%% === STEP 4 – Fit Group model
compareList = {'CTL','EPD','DYS'};
results_interaction = struct;

% === Load clinical data
fileClin = '...\clinical_data.xlsx';
Tclin = readtable(fileClin, 'VariableNamingRule', 'preserve');
Tclin.ID = replace(lower(string(Tclin.ID)), 'ta', 'x');  % normalize IDs

%% === STEP 4 – Robust fit “Group + Age + Sex” (morning only)
compareList = {'EPD','ADV','DYS'}; %to change
result_morning = struct;

% --- 4.0: string conversion and morning subset
Group_str   = string(Group);
Session_str = string(Session);
morningIdx  = (Session_str == referenceSession);
Group_m     = Group_str(morningIdx);
psd_m       = psd_theta(:, morningIdx);
Id_m        = Id(morningIdx);

% --- 4.1: extract clinical covariates for morning subjects
id_clean_m = lower(replace(string(Id_m), 'ta', 'x'));
[found_m, idxC_m] = ismember(id_clean_m, Tclin.ID);

age_m  = NaN(size(Id_m));
sex_m  = NaN(size(Id_m));
age_m(found_m)  = Tclin.Age(idxC_m(found_m));
sex_m(found_m)  = double(strcmpi(Tclin.Sex(idxC_m(found_m)),'M'));

nCh = size(psd_m,1);

for g = 1:numel(compareList)
    grp    = compareList{g};
    sel    = (Group_m == referenceGroup) | (Group_m == grp);
    nSubj  = sum(sel);
    fprintf('\n>>> Comparison %s vs %s (morning): %d subjects\n', ...
            referenceGroup, grp, nSubj);
    if nSubj < 4
        warning('  Too few observations, skipping %s vs %s.', referenceGroup, grp);
        continue;
    end
    
    % extract data
    y_mat    = psd_m(:, sel);           % [nCh x nSubj]
    Gcmp     = Group_m(sel);
    age_cmp  = age_m(sel);
    sex_cmp  = sex_m(sel);
    
    % imputation: median for age, mode for sex
    age_cmp  = fillmissing(age_cmp, 'constant', median(age_cmp, 'omitnan'));
    tmpSex   = sex_cmp(~isnan(sex_cmp));
    if ~isempty(tmpSex)
        sex_cmp = fillmissing(sex_cmp,'constant', mode(tmpSex));
    else
        sex_cmp = zeros(size(sex_cmp));
    end
    
    % preallocation
    tvals  = nan(nCh,1);
    pvals  = nan(nCh,1);
    d_vals = nan(nCh,1);
    
    for ch = 1:nCh
        y = y_mat(ch,:)';             % [nSubj x 1]
        G = double(Gcmp == grp);      % 0=CTL,1=grp
        
        % --- robust fit on [G, age, sex]
        X = [G, age_cmp, sex_cmp];    % [nSubj x 3]
        mdl = fitlm(X, y, 'RobustOpts','bisquare');
        
        CT = mdl.Coefficients;
        if height(CT) >= 2
            tvals(ch) = CT.tStat(2);
            pvals(ch) = CT.pValue(2);
        end
        
        % --- Cohen's d (pooled)
        y0 = y(G==0); y1 = y(G==1);
        if numel(y0)>1 && numel(y1)>1
            s0 = var(y0,1); s1 = var(y1,1);
            sp = sqrt(((numel(y0)-1)*s0 + (numel(y1)-1)*s1)/(numel(y0)+numel(y1)-2));
            d_vals(ch) = (mean(y1)-mean(y0))/sp;
        end
    end
    
    % save results
    result_morning.(grp).tvals    = tvals;
    result_morning.(grp).pvals    = pvals;
    result_morning.(grp).cohen_d  = d_vals;
    result_morning.(grp).chanlocs = chanlocs;
    result_morning.(grp).Id       = Id_m(sel);
end

%% === STEP 5 – Permutations + TFCE (Group+Age+Sex) ===
nPerms = 2000;
E_H    = [0.66 2];
progressStep = 100;

% Disable robustfit warning
warning('off','stats:statrobustfit:IterationLimit');

tfce_results = struct();
fprintf(' TFCE Permutations (%d permutations)\n\n', nPerms);

for g = 1:numel(compareList)
    grp = compareList{g};
    tag = sprintf('TFCE_%s_vs_%s', referenceGroup, grp);
    
    % subset morning data
    sel    = (Group_m == referenceGroup) | (Group_m == grp);
    y_mat  = psd_m(:, sel);
    Gcmp   = Group_m(sel);
    age_c  = age_m(sel);
    sex_c  = sex_m(sel);
    nCh    = size(y_mat,1);
    nSubj  = sum(sel);
    
    % imputation
    age_c = fillmissing(age_c, 'constant', median(age_c, 'omitnan'));
    tmpSex = sex_c(~isnan(sex_c));
    if ~isempty(tmpSex)
        sex_c = fillmissing(sex_c, 'constant', mode(tmpSex));
    else
        sex_c = zeros(size(sex_c));
    end
    
    % observed t-values
    t_obs = nan(nCh,1);
    for ch = 1:nCh
        y = y_mat(ch,:)';
        X = [double(Gcmp==grp), age_c, sex_c];
        mdl = fitlm(X, y, 'RobustOpts','bisquare');
        CT  = mdl.Coefficients;
        if height(CT)>=2, t_obs(ch) = CT.tStat(2); end
    end
    
    % permutations
    fprintf(' %d/%d completed...\n', nPerms, nPerms);
    t_perms = nan(nPerms, nCh);
    for p = 1:nPerms
        permIdx = randperm(nSubj);
        Gp      = double(Gcmp(permIdx)==grp);
        t_tmp   = nan(nCh,1);
        for ch = 1:nCh
            y   = y_mat(ch,:)';
            Xp  = [Gp, age_c, sex_c];
            mdlp = fitlm(Xp, y, 'RobustOpts','bisquare');
            CTp  = mdlp.Coefficients;
            if height(CTp)>=2, t_tmp(ch) = CTp.tStat(2); end
        end
        t_perms(p,:) = t_tmp';
        if mod(p,progressStep)==0 && p<nPerms
            rem = nPerms-p;
            fprintf(' %d/%d completed...\n', rem, nPerms);
        end
    end
    fprintf('   Permutations completed.\n');
    
    % TFCE computation
    Results = lmeEEG_TFCE_mod(t_obs, t_perms, chanlocs, E_H);
    nSig    = sum(Results.Mask);
    
    % summary
    fprintf(' T_Obs: [%d 1], T_Perms: [%d %d]\n', numel(t_obs), size(Results.TFCE_Perms));
    fprintf('Max TFCE=%.2f | Significant channels= %d\n\n', Results.maxTFCE, nSig);
    
    tfce_results.(tag).T_Obs      = t_obs;
    tfce_results.(tag).TFCE_Obs   = Results.TFCE_Obs;
    tfce_results.(tag).P_Values   = Results.P_Values;
    tfce_results.(tag).Mask       = Results.Mask;
    tfce_results.(tag).TFCE_Perms = Results.TFCE_Perms;
    tfce_results.(tag).maxTFCE    = Results.maxTFCE;
end

% Restore warnings
warning('on','stats:statrobustfit:IterationLimit');

%% === FIGURE 2×3 WITH FIXED COLOR SCALES FOR p-values ===
% Pre-compute common TFCE limit
allTF = [];
for g = 1:numel(compareList)
    grp = compareList{g};
    tag = sprintf('TFCE_%s_vs_%s', referenceGroup, grp);
    allTF = [allTF; tfce_results.(tag).TFCE_Obs];
end
mxTF = max(abs(allTF));

% Create figure
figure('Color','w','Name','TFCE – Main Effect Group','Units','normalized','Position',[.1 .1 .8 .6]);
t = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
title(t, sprintf('TFCE – Main Effect Group  |  Baseline: %s – %s', referenceGroup, referenceSession), ...
      'FontSize',14,'FontWeight','bold');

for col = 1:3
    grp = compareList{col};
    tag = sprintf('TFCE_%s_vs_%s', referenceGroup, grp);
    res = tfce_results.(tag);
    
    % TFCE map
    nexttile(col);
    topoplot(res.TFCE_Obs, chanlocs, ...
        'maplimits',[-mxTF mxTF], ...
        'electrodes','on','style','both');
    colorbar;
    title(sprintf('%s vs %s — TFCE', grp, referenceGroup), 'FontSize',10);
    
    % –log10(p) map with fixed scale [0 3]
    nexttile(col+3);
    logp  = -log10(res.P_Values);
    sigCh = find(res.Mask);
    topoplot(logp, chanlocs, ...
        'maplimits',[0 3], ...
        'electrodes','on','style','both', ...
        'emarker2',{sigCh,'o',[1 0 0],4,2});
    colorbar;
    title(sprintf('%s vs %s — -log_{10}(p)', grp, referenceGroup), 'FontSize',10);
end
