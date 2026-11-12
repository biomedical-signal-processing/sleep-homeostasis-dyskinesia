%% ========================================================================
% DYS — MINIMAL DIFF CORRELATIONS
%
% Analyses:
%   1) Delta SampEn  (Late - Early) vs clinical variable
%   2) Delta SWA     (Late - Early) vs clinical variable
%   3) Partial correlation: Delta SampEn vs clinical variable, controlling
%      for Delta SWA
%   4) EEG–EEG correlation: Delta SampEn vs Delta SWA
%   5) Optional regression: clinical variable ~ Delta SampEn + Delta SWA
%      (+ Age + Sex)
%
% Inputs:
%   - Subject-level SampEn matrices
%   - Subject-level scoring files
%   - SWA globalRes output
%   - Clinical spreadsheet
%
% Outputs:
%   - 3-panel clinical correlation figure
%   - EEG–EEG figure
%   - Correlation summary CSV
%   - Optional regression summary TXT + CSV
% ========================================================================

clear; clc; close all; rng(42);

%% ===================== PATHS ============================================
SWTEST_MFILE = '...\swtest\swtest.m';
addpath(fileparts(SWTEST_MFILE));

ROOT        = '...\Annals_processed';
SCORING_DIR = fullfile(ROOT, 'infoscoringpersubject');
MAP_PATH    = '...\map_149to128.mat';

CLIN_XLSX   = '...\db_annals_theta_CLEAN.xlsx';
CLIN_SHEET  = 1;

SWA_GLOBALRES = '...\globalRes_SWA_CLAMP.mat';

OUT_DIR = '...\GroupAnalysis\DYS_corr_MINIMAL_DIFF';
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

%% ===================== USER OPTIONS =====================================
% SampEn input
METRIC_KIND    = 'SampleEntropy';
METRIC_SUBDIR  = 'SampEn_output';
METRIC_VARNAME = '';      % empty = autodetect
MSE_SCALE      = 1;

% Sleep state selection
STATE     = 'N3';         % 'N3' or 'N2N3'
NSEG      = 10;
EARLY_SEG = 2:5;
LATE_SEG  = 6:9;

% Epoch aggregation
AGG_METHOD   = 'mean';    % 'mean' or 'median'
LOG_STRATEGY = 'post';    % 'post', 'pre', or 'none'

% Epoch balancing
EPOCH_BALANCING    = 'equal_per_subject';
EPOCH_SELECT_MODE  = 'first';
MIN_EPOCHS_PER_SEG = 1;
RNG_SEED_BALANCING = 1234;
rng(RNG_SEED_BALANCING);

% Clinical variable
CLIN_VAR  = 'UPDRS_IV_II';
CLIN_TITL = 'MDS-UPDRS IV item 2';

% Spatial selection
SPATIAL_MODE = 'roi';    % 'all' | 'roi' | 'channels' | 'single'
ROI_NAME     = 'f';
CHAN_LIST    = {};

% SWA input
SWA_GROUP   = 'DYS';
SWA_BAND    = 'B_0p5_4';
SWA_FEATURE = 'event';
LOG_SWA     = 'log';          % 'log' or 'none'

% Optional regression
DO_REGRESSION   = true;
INCLUDE_AGE_SEX = false;

% Figure settings
FIG_W = 1200;
FIG_H = 420;
DPI   = 180;

fprintf('\n================= START DYS MINIMAL DIFF CORRELATIONS =================\n');
fprintf('SampEn: %s | state=%s | log=%s | agg=%s\n', METRIC_SUBDIR, STATE, LOG_STRATEGY, AGG_METHOD);
fprintf('SWA: %s | %s.%s\n', SWA_GLOBALRES, SWA_BAND, SWA_FEATURE);
fprintf('Clinical variable: %s (%s)\n', CLIN_VAR, CLIN_TITL);
fprintf('Output: %s\n\n', OUT_DIR);

%% ===================== LOAD CHANNEL MAPPING ==============================
keep_idx_149to128 = [];
kept128_labels    = strings(0,1);

if exist(MAP_PATH, 'file') == 2
    Smap = load(MAP_PATH);
    if isfield(Smap, 'keep_idx_149to128')
        keep_idx_149to128 = Smap.keep_idx_149to128(:);
    end
    if isfield(Smap, 'kept128_labels')
        kept128_labels = upper(string(Smap.kept128_labels(:)'));
    end
    fprintf('[MAP] Loaded mapping: %d indices, %d labels\n', numel(keep_idx_149to128), numel(kept128_labels));
else
    fprintf('[MAP] Mapping file not found. Using matrix as-is.\n');
end

%% ===================== LOAD CLINICAL TABLE ==============================
Tdb = readtable(CLIN_XLSX, 'Sheet', CLIN_SHEET, 'VariableNamingRule', 'preserve');
Tdb.ID    = upper(string(Tdb.ID));
Tdb.Group = lower(string(Tdb.Group));
Tdb.Sex   = upper(string(Tdb.Sex));
Tdb.Sex01 = sex_to01(Tdb.Sex);

requiredCols = {'ID','Group','Age','Sex','Sex01',CLIN_VAR};
missingCols  = setdiff(requiredCols, Tdb.Properties.VariableNames);
if ~isempty(missingCols)
    error('Clinical file is missing columns: %s', strjoin(missingCols, ', '));
end

Tdb = Tdb(strcmpi(Tdb.Group, 'dys'), :);
if isempty(Tdb)
    error('No DYS subjects found in the clinical file.');
end

fprintf('[CLIN] DYS rows in clinical table: %d\n', height(Tdb));

%% ===================== BUILD SAMPEN TABLE ===============================
Tmet = table('Size', [0 7], ...
    'VariableTypes', {'string','string','double','double','double','double','double'}, ...
    'VariableNames', {'ID','Group','EARLY','LATE','DIFF','Age','Sex01'});

listing = get_subject_listing(ROOT, METRIC_SUBDIR, 'DYS');
fprintf('[SampEn] Files found: %d\n', numel(listing));

nSkip = 0;

for k = 1:numel(listing)
    subj    = listing(k).subj;
    matPath = listing(k).fullpath;

    try
        [M2D, labels_in] = load_metric2D_with_labels(matPath, METRIC_VARNAME, METRIC_KIND, MSE_SCALE);
        [M2D, labels]    = map_or_label(M2D, labels_in, keep_idx_149to128, kept128_labels);

        idxSel = select_spatial_indices(SPATIAL_MODE, ROI_NAME, CHAN_LIST, labels);
        if isempty(idxSel)
            idxSel = 1:size(M2D,1);
        end
        Msel = M2D(idxSel, :);

        [epochStages, okScore] = load_scoring(SCORING_DIR, 'DYS', subj, size(Msel,2));
        if ~okScore
            nSkip = nSkip + 1;
            continue;
        end

        [Esel, Lsel, okSel] = select_E_L_indices_strict(epochStages, NSEG, EARLY_SEG, LATE_SEG, STATE);
        if ~okSel
            nSkip = nSkip + 1;
            continue;
        end

        [Esel, Lsel, okClamp] = apply_epoch_clamp( ...
            Esel, Lsel, EPOCH_BALANCING, EPOCH_SELECT_MODE, MIN_EPOCHS_PER_SEG, RNG_SEED_BALANCING);
        if ~okClamp
            nSkip = nSkip + 1;
            continue;
        end

        % Collapse selected channels to one epoch-wise vector
        E_ep = squeeze(mean(Msel(:, Esel), 1, 'omitnan'));
        L_ep = squeeze(mean(Msel(:, Lsel), 1, 'omitnan'));

        % Aggregate and optionally log-transform
        switch lower(LOG_STRATEGY)
            case 'pre'
                E_ep = safe_log(E_ep);
                L_ep = safe_log(L_ep);
                EAR = agg(E_ep, AGG_METHOD);
                LAT = agg(L_ep, AGG_METHOD);

            case 'post'
                EAR = agg(E_ep, AGG_METHOD);
                LAT = agg(L_ep, AGG_METHOD);
                EAR = safe_log(EAR);
                LAT = safe_log(LAT);

            case 'none'
                EAR = agg(E_ep, AGG_METHOD);
                LAT = agg(L_ep, AGG_METHOD);

            otherwise
                error('LOG_STRATEGY must be post, pre, or none.');
        end

        DIFF = LAT - EAR;

        idxC = find(Tdb.ID == upper(string(subj)), 1, 'first');
        if isempty(idxC)
            nSkip = nSkip + 1;
            continue;
        end

        Age   = double(Tdb.Age(idxC));
        Sex01 = double(Tdb.Sex01(idxC));

        Tmet = [Tmet; {upper(string(subj)),'dys',EAR,LAT,DIFF,Age,Sex01}]; %#ok<AGROW>

    catch ME
        nSkip = nSkip + 1;
        warning('Skipping %s: %s', subj, ME.message);
    end
end

if isempty(Tmet)
    error('No valid SampEn subjects after preprocessing.');
end

fprintf('[SampEn] Final table: N=%d (skipped=%d)\n', height(Tmet), nSkip);

%% ===================== BUILD SWA TABLE ==================================
Tswa = build_swa_table_globalRes( ...
    SWA_GLOBALRES, SWA_GROUP, SWA_BAND, SWA_FEATURE, ...
    SPATIAL_MODE, ROI_NAME, CHAN_LIST, LOG_SWA);

fprintf('[SWA] Final table: N=%d\n', height(Tswa));

%% ===================== MERGE TABLES =====================================
Tclin = Tdb(:, {'ID','Age','Sex01',CLIN_VAR});

Tjoin = innerjoin(Tclin, Tmet(:, {'ID','DIFF'}), 'Keys', 'ID');
Tjoin = innerjoin(Tjoin, Tswa(:, {'ID','SWA_DIFF'}), 'Keys', 'ID');

fprintf('[JOIN] Final N=%d\n', height(Tjoin));
disp('[JOIN IDs]');
disp(Tjoin.ID(:)');

if height(Tjoin) < 3
    error('N < 3 after merge. Check ID matching across sources.');
end

%% ===================== CORRELATIONS =====================================
x     = double(Tjoin.(CLIN_VAR));
dSamp = double(Tjoin.DIFF);
dSWA  = double(Tjoin.SWA_DIFF);

% 1) Delta SampEn vs clinical
[method1, r1, p1, N1, hNx1, pNx1, hNy1, pNy1] = corr_auto_swtest(x, dSamp);

% 2) Delta SWA vs clinical
[method2, r2, p2, N2, hNx2, pNx2, hNy2, pNy2] = corr_auto_swtest(x, dSWA);

% 3) Partial: Delta SampEn vs clinical controlling for Delta SWA
[xr, yr] = residualize_pair_for_plot(x, dSamp, dSWA);
[method3, r3, p3, N3, hNxr, pNxr, hNyr, pNyr] = partialcorr_auto_swtest(x, dSamp, dSWA);

% 4) EEG–EEG: Delta SampEn vs Delta SWA
[methodEEG, rEEG, pEEG, NEEG, hNxE, pNxE, hNyE, pNyE] = corr_auto_swtest(dSamp, dSWA);

fprintf('\n[DIFF] Delta SampEn vs %s | %s r=%.3f p=%.3g N=%d\n', CLIN_VAR, method1, r1, p1, N1);
fprintf('[DIFF] Delta SWA   vs %s | %s r=%.3f p=%.3g N=%d\n', CLIN_VAR, method2, r2, p2, N2);
fprintf('[DIFF] Partial controlling Delta SWA | %s r=%.3f p=%.3g N=%d\n', method3, r3, p3, N3);
fprintf('[EEG–EEG] Delta SampEn vs Delta SWA | %s r=%.3f p=%.3g N=%d\n', methodEEG, rEEG, pEEG, NEEG);

%% ===================== FIGURE: 3 CLINICAL PANELS ========================
fig = figure('Color', 'w', 'Position', [120 120 FIG_W FIG_H]);
tl  = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, 'PDdys^{+}: clinical correlations on Late–Early differences', 'FontWeight', 'bold');

ax1 = nexttile(tl, 1);
plot_scatter_fit(ax1, x, dSamp, CLIN_TITL, 'Delta SampEn (Late - Early)', ...
    sprintf('%s r=%.2f, p=%.3g, N=%d', method1, r1, p1, N1));

ax2 = nexttile(tl, 2);
plot_scatter_fit(ax2, x, dSWA, CLIN_TITL, 'Delta SWA (Late - Early)', ...
    sprintf('%s r=%.2f, p=%.3g, N=%d', method2, r2, p2, N2));

ax3 = nexttile(tl, 3);
plot_scatter_fit(ax3, xr, yr, ...
    sprintf('%s residualized for Delta SWA', CLIN_TITL), ...
    'Delta SampEn residualized for Delta SWA', ...
    sprintf('partialcorr (%s) r=%.2f, p=%.3g, N=%d', method3, r3, p3, N3));

figClin = fullfile(OUT_DIR, sprintf('DYS_DIFF_corr_3panels_vs_%s.png', lower(CLIN_VAR)));
exportgraphics(fig, figClin, 'Resolution', DPI);

%% ===================== FIGURE: EEG–EEG ==================================
figEEG = figure('Color', 'w', 'Position', [120 120 round(FIG_W*0.55) FIG_H]);
axEEG = axes(figEEG); %#ok<LAXES>
plot_scatter_fit(axEEG, dSamp, dSWA, ...
    'Delta SampEn (Late - Early)', ...
    'Delta SWA (Late - Early)', ...
    sprintf('EEG–EEG: %s r=%.2f, p=%.3g, N=%d', methodEEG, rEEG, pEEG, NEEG));

figEEGFile = fullfile(OUT_DIR, 'DYS_DIFF_EEG_EEG_dSampEn_vs_dSWA.png');
exportgraphics(figEEG, figEEGFile, 'Resolution', DPI);

%% ===================== OPTIONAL REGRESSION ==============================
regTable = table();
if DO_REGRESSION
    ok = isfinite(x) & isfinite(dSamp) & isfinite(dSWA);

    if INCLUDE_AGE_SEX
        Age0 = double(Tjoin.Age);
        Sex0 = double(Tjoin.Sex01);
        ok = ok & isfinite(Age0) & isfinite(Sex0);
    end

    y0    = x(ok);
    samp0 = dSamp(ok);
    swa0  = dSWA(ok);

    if INCLUDE_AGE_SEX
        Age0 = double(Tjoin.Age(ok));
        Sex0 = double(Tjoin.Sex01(ok));
        Treg = table(y0, samp0, swa0, Age0, Sex0, ...
            'VariableNames', {'Clinical','dSampEn','dSWA','Age','Sex01'});
        mdl = fitlm(Treg, 'Clinical ~ dSampEn + dSWA + Age + Sex01');
        VIF = calc_vif([samp0 swa0 Age0 Sex0]);
    else
        Treg = table(y0, samp0, swa0, ...
            'VariableNames', {'Clinical','dSampEn','dSWA'});
        mdl = fitlm(Treg, 'Clinical ~ dSampEn + dSWA');
        VIF = calc_vif([samp0 swa0]);
    end

    C = mdl.Coefficients;

    beta_samp = C.Estimate(strcmp(C.Properties.RowNames, 'dSampEn'));
    p_samp    = C.pValue (strcmp(C.Properties.RowNames, 'dSampEn'));
    beta_swa  = C.Estimate(strcmp(C.Properties.RowNames, 'dSWA'));
    p_swa     = C.pValue (strcmp(C.Properties.RowNames, 'dSWA'));

    fprintf('\n[REG] Clinical ~ dSampEn + dSWA%s | N=%d\n', ...
        ternary(INCLUDE_AGE_SEX, ' + Age + Sex', ''), numel(y0));
    fprintf('[REG] beta(dSampEn)=%.4g p=%.3g | beta(dSWA)=%.4g p=%.3g\n', ...
        beta_samp, p_samp, beta_swa, p_swa);
    fprintf('[REG] VIF: %s\n', mat2str(VIF, 3));

    txt = evalc('disp(mdl)');
    regTxt = fullfile(OUT_DIR, sprintf('DYS_DIFF_regression_%s.txt', ternary(INCLUDE_AGE_SEX, 'AgeSex', 'basic')));
    fid = fopen(regTxt, 'w');
    fprintf(fid, '%s', txt);
    fclose(fid);

    regTable = table( ...
        numel(y0), INCLUDE_AGE_SEX, beta_samp, p_samp, beta_swa, p_swa, ...
        {mat2str(VIF,3)}, ...
        'VariableNames', {'N','IncludeAgeSex','Beta_dSampEn','P_dSampEn','Beta_dSWA','P_dSWA','VIF'});
end

%% ===================== SAVE CORRELATION SUMMARY ==========================
CorrSummary = table();

CorrSummary.Analysis = [ ...
    "Delta SampEn vs Clinical"; ...
    "Delta SWA vs Clinical"; ...
    "Partial Delta SampEn vs Clinical | control Delta SWA"; ...
    "EEG–EEG Delta SampEn vs Delta SWA"];

CorrSummary.Method = [string(method1); string(method2); string(method3); string(methodEEG)];
CorrSummary.R      = [r1; r2; r3; rEEG];
CorrSummary.P      = [p1; p2; p3; pEEG];
CorrSummary.N      = [N1; N2; N3; NEEG];

CorrSummary.swtest_x_h = [hNx1; hNx2; hNxr; hNxE];
CorrSummary.swtest_x_p = [pNx1; pNx2; pNxr; pNxE];
CorrSummary.swtest_y_h = [hNy1; hNy2; hNyr; hNyE];
CorrSummary.swtest_y_p = [pNy1; pNy2; pNyr; pNyE];

CorrSummary.ClinicalVar = repmat(string(CLIN_VAR), height(CorrSummary), 1);
CorrSummary.ClinicalLabel = repmat(string(CLIN_TITL), height(CorrSummary), 1);

corrCSV = fullfile(OUT_DIR, 'DYS_DIFF_corr_summary.csv');
writetable(CorrSummary, corrCSV);

if ~isempty(regTable)
    regCSV = fullfile(OUT_DIR, 'DYS_DIFF_regression_summary.csv');
    writetable(regTable, regCSV);
end

fprintf('\n✔ Clinical figure: %s\n', figClin);
fprintf('✔ EEG–EEG figure: %s\n', figEEGFile);
fprintf('✔ Correlation CSV: %s\n', corrCSV);
if ~isempty(regTable)
    fprintf('✔ Regression CSV: %s\n', regCSV);
end
fprintf('================== DONE ==================\n');

%% ========================================================================
%%                              LOCAL FUNCTIONS
%% ========================================================================

function Tswa = build_swa_table_globalRes(matPath, groupName, bandName, featName, spatialMode, roiName, chanList, logMode)
S = load(matPath);
GR = S.globalRes;

if isfield(GR, groupName)
    G = GR.(groupName);
elseif isfield(GR, lower(groupName))
    G = GR.(lower(groupName));
else
    error('Group "%s" not found in SWA globalRes.', groupName);
end

IDs    = upper(string(G.subj(:)));
labels = upper(string(G.labels(:)));

idxSel = select_spatial_indices(spatialMode, roiName, chanList, labels);

E = force_numeric_mat(G.(bandName).early.(featName), 'SWA early');
L = force_numeric_mat(G.(bandName).late.(featName),  'SWA late');

if size(E,1) ~= numel(labels) || size(L,1) ~= numel(labels)
    error('SWA matrix size does not match SWA labels.');
end
if size(E,2) ~= numel(IDs) || size(L,2) ~= numel(IDs)
    error('SWA subject count mismatch.');
end

EAR = mean(E(idxSel,:), 1, 'omitnan')';
LAT = mean(L(idxSel,:), 1, 'omitnan')';

if strcmpi(logMode, 'log')
    EAR = safe_log(EAR);
    LAT = safe_log(LAT);
end

DIFF = LAT - EAR;
Tswa = table(IDs, EAR, LAT, DIFF, 'VariableNames', {'ID','SWA_EARLY','SWA_LATE','SWA_DIFF'});
end

function X = force_numeric_mat(V, tag)
if isnumeric(V) && ismatrix(V)
    X = double(V);
elseif isstruct(V) && isfield(V, 'mat') && isnumeric(V.mat)
    X = double(V.mat);
else
    error('%s is not a numeric matrix.', tag);
end
end

function listing = get_subject_listing(ROOT, metricSubdir, G)
p = fullfile(ROOT, metricSubdir, G);
if ~exist(p, 'dir')
    listing = [];
    return;
end

d = dir(fullfile(p, '*_eeg_minimal_*.mat'));
listing = struct('subj',{},'fullpath',{},'name',{});

for i = 1:numel(d)
    subj = extract_subj_id(d(i).name);
    listing(end+1).subj   = subj; %#ok<AGROW>
    listing(end).fullpath = fullfile(d(i).folder, d(i).name);
    listing(end).name     = d(i).name;
end
end

function subj = extract_subj_id(fname)
tok = regexp(fname, '^(.*?)_eeg_minimal', 'tokens', 'once');
if isempty(tok)
    subj = erase(fname, '.mat');
else
    subj = tok{1};
end
end

function [M2D, labels] = load_metric2D_with_labels(matPath, metricVarName, metricKind, MSE_SCALE)
S = load(matPath);
labels = get_labels_from_struct(S);
M2D = [];

    function X = pick2D(Xin, scale)
        if ndims(Xin) == 2
            X = Xin;
        elseif ndims(Xin) == 3
            [~,~,K] = size(Xin);
            assert(scale >= 1 && scale <= K, 'MSE_SCALE out of range.');
            X = Xin(:,:,scale);
        else
            error('Unsupported number of dimensions.');
        end
    end

if ~isempty(metricVarName) && isfield(S, metricVarName) && isnumeric(S.(metricVarName))
    M2D = pick2D(S.(metricVarName), MSE_SCALE);
else
    switch lower(string(metricKind))
        case "sampleentropy"
            cand = {'SampleEntropy','SampEn'};
        case "spectralentropy"
            cand = {'SpectralEntropy','SpecEn'};
        case "msesampen"
            cand = {'MSESampEn','MSE_SampEn'};
        case "msespecen"
            cand = {'MSESpecEn','MSE_SpecEn','MSE_SpectralEntropy'};
        otherwise
            cand = {};
    end

    for i = 1:numel(cand)
        if isfield(S, cand{i}) && isnumeric(S.(cand{i}))
            M2D = pick2D(S.(cand{i}), MSE_SCALE);
            break;
        end
    end

    if isempty(M2D)
        fns = fieldnames(S);
        for i = 1:numel(fns)
            V = S.(fns{i});
            if isnumeric(V) && (ismatrix(V) || ndims(V)==3)
                M2D = pick2D(V, MSE_SCALE);
                break;
            end
        end
    end
end

if isempty(M2D)
    error('Metric not found in MAT file.');
end

if isempty(labels)
    labels = "E" + (1:size(M2D,1));
end
labels = upper(string(labels));
end

function labels = get_labels_from_struct(S)
labels = string([]);
cand = {'kept128_labels','chanlabels','labels','ChLabels','Channels'};
for i = 1:numel(cand)
    if isfield(S, cand{i})
        labels = string(S.(cand{i})(:).');
        break;
    end
end
labels = upper(string(labels));
end

function [Mout, labelsOut] = map_or_label(Min, labels_in, keep_idx_149to128, kept128_labels)
nCh = size(Min,1);
labelsOut = labels_in;

if nCh == 149 && ~isempty(keep_idx_149to128)
    Mout = Min(keep_idx_149to128, :);
    if ~isempty(kept128_labels)
        labelsOut = kept128_labels;
    end
else
    Mout = Min;
    if isempty(labelsOut)
        labelsOut = "E" + (1:nCh);
    end
end

labelsOut = upper(string(labelsOut));
end

function idxSel = select_spatial_indices(mode, ROI_NAME, CHAN_LIST, labels)
switch lower(mode)
    case 'all'
        idxSel = 1:numel(labels);

    case 'roi'
        ROI = roi_labels();
        rn  = lower(string(ROI_NAME));
        assert(isfield(ROI, rn), 'ROI "%s" not recognized.', ROI_NAME);
        wanted = upper(string(ROI.(rn)));
        idxSel = labels_match_indices(wanted, labels);

    case 'channels'
        assert(~isempty(CHAN_LIST), 'CHAN_LIST is empty.');
        wanted = upper(string(CHAN_LIST));
        idxSel = labels_match_indices(wanted, labels);

    case 'single'
        assert(~isempty(CHAN_LIST), 'CHAN_LIST is empty.');
        wanted = upper(string(CHAN_LIST(1)));
        idxSel = labels_match_indices(wanted, labels);
        if numel(idxSel) > 1
            idxSel = idxSel(1);
        end

    otherwise
        error('Invalid SPATIAL_MODE: %s', mode);
end
end

function idx = labels_match_indices(wanted, labels)
L = upper(string(labels));
idx = [];

for i = 1:numel(wanted)
    w = upper(string(wanted(i)));
    w = regexprep(w, '^E0*', 'E');
    j = find(L == w, 1, 'first');
    if ~isempty(j)
        idx(end+1) = j; %#ok<AGROW>
    end
end

idx = unique(idx);
end


function [epochStages, ok] = load_scoring(SCORING_DIR, G, subj, nEpoch)
ok = true;
scorePath = fullfile(SCORING_DIR, sprintf('%s_%s_scoring.mat', G, subj));

if exist(scorePath, 'file') ~= 2
    ok = false;
    epochStages = [];
    return;
end

T = load(scorePath);
epochStages = T.epochStages(:)';

if numel(epochStages) ~= nEpoch
    warning('Scoring mismatch for %s: %d vs %d', subj, numel(epochStages), nEpoch);
    ok = false;
end
end

function [Esel, Lsel, ok] = select_E_L_indices_strict(epochStages, NSEG, earlySeg, lateSeg, state_sel)
ok = false;
Esel = [];
Lsel = [];

idxNREM = find(epochStages==1 | epochStages==2 | epochStages==3);
if numel(idxNREM) < NSEG
    return;
end

edges = round(linspace(1, numel(idxNREM)+1, NSEG+1));

Eidx = [];
Lidx = [];
for s = earlySeg
    Eidx = [Eidx, idxNREM(edges(s):edges(s+1)-1)]; %#ok<AGROW>
end
for s = lateSeg
    Lidx = [Lidx, idxNREM(edges(s):edges(s+1)-1)]; %#ok<AGROW>
end

switch upper(state_sel)
    case 'N2N3'
        mask = (epochStages==2 | epochStages==3);
    case 'N3'
        mask = (epochStages==3);
    otherwise
        error('STATE must be N2N3 or N3.');
end

Esel = Eidx(mask(Eidx));
Lsel = Lidx(mask(Lidx));

if isempty(Esel) || isempty(Lsel)
    return;
end
ok = true;
end

function [Esel, Lsel, ok] = apply_epoch_clamp(Esel, Lsel, mode, pick, min_each, seed)
if nargin < 6, seed = 1234; end
ok = true;

if strcmpi(mode, 'none')
    if min(numel(Esel), numel(Lsel)) < min_each
        ok = false;
    end
    return;
end

m = min(numel(Esel), numel(Lsel));
if m < min_each
    ok = false;
    return;
end

switch lower(pick)
    case 'first'
        Esel = Esel(1:m);
        Lsel = Lsel(1:m);

    case 'last'
        Esel = Esel(end-m+1:end);
        Lsel = Lsel(end-m+1:end);

    case 'center'
        Esel = center_block(Esel, m);
        Lsel = center_block(Lsel, m);

    case 'random'
        rng(seed);
        Esel = Esel(randperm(numel(Esel), m));
        Lsel = Lsel(randperm(numel(Lsel), m));

    otherwise
        error('Invalid EPOCH_SELECT_MODE.');
end
end

function ids = center_block(pool, m)
N = numel(pool);
start = floor((N-m)/2) + 1;
ids = pool(start:start+m-1);
end

function y = safe_log(x)
x = double(x);
x(x <= 0) = eps;
y = log(x);
end

function out = agg(v, how)
if strcmpi(how, 'median')
    out = median(v, 'omitnan');
else
    out = mean(v, 'omitnan');
end
end

function s01 = sex_to01(SexStr)
s = upper(string(SexStr));
s01 = NaN(size(s));
s01(s=="M") = 1;
s01(s=="F") = 0;
end

function plot_scatter_fit(ax, x, y, xlab, ylab, ttl)
ok = isfinite(x) & isfinite(y);
x = double(x(ok));
y = double(y(ok));

if numel(x) < 3 || range(x)==0 || range(y)==0
    axis(ax, 'off');
    text(ax, 0.5, 0.5, 'N < 3 or zero variance', 'HorizontalAlignment', 'center');
    return;
end

scatter(ax, x, y, 60, 'filled');
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

xlabel(ax, xlab, 'FontWeight', 'bold');
ylabel(ax, ylab, 'FontWeight', 'bold');
title(ax, ttl, 'FontWeight', 'bold');

mdl = fitlm(x, y);
xg  = linspace(min(x), max(x), 200)';
[yh, yci] = predict(mdl, xg);

fill(ax, [xg; flipud(xg)], [yci(:,1); flipud(yci(:,2))], 0.9*[1 1 1], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);
plot(ax, xg, yh, 'LineWidth', 2.0);
end

function vif = calc_vif(X)
X = double(X);
k = size(X,2);
vif = NaN(1,k);

for j = 1:k
    y = X(:,j);
    Xo = X(:, setdiff(1:k, j));
    ok = all(isfinite([y Xo]), 2);
    if nnz(ok) < 3
        continue;
    end
    mdl = fitlm(Xo(ok,:), y(ok));
    R2 = mdl.Rsquared.Ordinary;
    vif(j) = 1 / max(eps, 1 - R2);
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function [method, r, p, N, hx, px, hy, py] = corr_auto_swtest(x, y)
x = double(x(:));
y = double(y(:));

ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
N = numel(x);

r = NaN; p = NaN; method = 'NA';
hx = NaN; px = NaN; hy = NaN; py = NaN;

if N < 3 || range(x)==0 || range(y)==0
    return;
end

[hx, px] = swtest_safe(x);
[hy, py] = swtest_safe(y);

if hx==0 && hy==0
    method = 'Pearson';
    [r, p] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
else
    method = 'Spearman';
    [r, p] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
end
end

function [method, r, p, N, hxr, pxr, hyr, pyr] = partialcorr_auto_swtest(x, y, z)
x = double(x(:));
y = double(y(:));
z = double(z(:));

ok = isfinite(x) & isfinite(y) & isfinite(z);
x = x(ok);
y = y(ok);
z = z(ok);
N = numel(x);

r = NaN; p = NaN; method = 'NA';
hxr = NaN; pxr = NaN; hyr = NaN; pyr = NaN;

if N < 3 || range(x)==0 || range(y)==0 || range(z)==0
    return;
end

xr = residualize_1cov(x, z);
yr = residualize_1cov(y, z);

[hxr, pxr] = swtest_safe(xr);
[hyr, pyr] = swtest_safe(yr);

if hxr==0 && hyr==0
    method = 'Pearson';
    [r, p] = partialcorr(x, y, z, 'Type', 'Pearson', 'Rows', 'complete');
else
    method = 'Spearman';
    [r, p] = partialcorr(x, y, z, 'Type', 'Spearman', 'Rows', 'complete');
end
end

function [xr, yr] = residualize_pair_for_plot(x, y, z)
xr = residualize_1cov(x, z);
yr = residualize_1cov(y, z);
end

function [h, p] = swtest_safe(v)
v = double(v(:));
v = v(isfinite(v));

if numel(v) < 3
    h = 1;
    p = NaN;
    return;
end

try
    [h, p] = swtest(v);
catch
    h = 1;
    p = NaN;
end
end

function r = residualize_1cov(y, c)
y = double(y(:));
c = double(c(:));

ok = isfinite(y) & isfinite(c);
r = NaN(size(y));

if nnz(ok) < 3 || range(c(ok)) == 0
    r(ok) = y(ok) - mean(y(ok), 'omitnan');
    return;
end

X = [ones(nnz(ok),1) c(ok)];
b = X \ y(ok);
r(ok) = y(ok) - X*b;
end