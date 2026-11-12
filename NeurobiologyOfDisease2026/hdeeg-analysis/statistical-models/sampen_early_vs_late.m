%% ========================================================================
% SAMPEN — EARLY vs LATE (WITHIN-GROUP) — 128 CHANNELS — TFCE
%
% This script compares Early vs Late sleep within each group using
% sample entropy values already computed at the epoch level.
%
% Pipeline:
%   1) Load subject-wise SampEn matrices
%   2) Load epoch staging
%   3) Select Early and Late epochs from NREM-based segmentation
%   4) Keep only the requested state (e.g., N3)
%   5) Clamp Early/Late to the same number of epochs per subject
%   6) Aggregate epochs to one value per subject and channel
%   7) Run paired t-statistics + sign-flip permutations + TFCE
%   8) Save CSV results and summary figures
%
% Requirements:
%   - lmeEEG_TFCE_mod on path
%   - EEGLAB topoplot for topographic plots
% ========================================================================

clear; clc; close all; rng(42);

%% ====================== SETTINGS ========================================
GROUPS = {'CTL','ADV','DYS'};
ANALYSIS_STATES = {'N3'};              % {'N3'} or {'N2N3','N3'}

% NREM segmentation
NSEG      = 10;
EARLY_SEG = 2:5;
LATE_SEG  = 6:9;

% Epoch aggregation
AGG_METHOD = 'mean';                   % 'mean' or 'median'

% Log transform
LOG_STRATEGY = 'post';                 % 'post', 'pre', or 'none'
EPS_LOG      = 1e-8;

% Epoch balancing
EPOCH_BALANCING    = 'equal_per_subject';   % 'equal_per_subject' or 'none'
EPOCH_SELECT_MODE  = 'first';               % 'first', 'last', 'center', 'random'
MIN_EPOCHS_PER_SEG = 3;
RNG_SEED_BALANCING = 1234;

% Statistics
TFCE_EH       = [0.66 2];
N_PERM        = 1000;
SIGNFLIP_MODE = 'per_channel';         % 'per_channel' or 'per_subject_coherent'

% Input / output
ROOT        = '...\_processed';
SCORING_DIR = fullfile(ROOT, 'infoscoringpersubject');
OUTROOT     = '...\_entropy';
OUT_DIR     = fullfile(OUTROOT, 'GroupAnalysis', 'WITHIN_SampEn');
ensure_dir(OUT_DIR);

% SampEn files
METRIC_SUBDIR = 'SampEn_output';
METRIC_NAME   = '';                    % '' = auto-detect SampleEntropy / SampEn

% Channel mapping: 149 -> 128
MAP_PATH   = '...\map_149to128.mat';   % keep_idx_149to128, kept128_labels
CHANLOCS_M = '...\chanlocs128.mat';    % chanlocs128

% Plotting
DPI = 300;
PRINT_COUNTS = true;

%% ====================== LOAD STATIC DATA ================================
Smap = load(MAP_PATH);
keep_idx_149to128 = Smap.keep_idx_149to128(:);
kept128_labels    = string(Smap.kept128_labels(:)');
L = load(CHANLOCS_M);
chanlocs = L.chanlocs128;

if exist('lmeEEG_TFCE_mod','file') ~= 2
    error('lmeEEG_TFCE_mod not found on MATLAB path.');
end
if exist('topoplot','file') ~= 2
    warning('topoplot not found. Topographic plots will use scatter fallback.');
end

fprintf('\n=== WITHIN-GROUP SAMPEN ===\n');
fprintf('LOG=%s | AGG=%s | CLAMP=%s/%s | N_PERM=%d\n', ...
    upper(LOG_STRATEGY), upper(AGG_METHOD), ...
    upper(EPOCH_BALANCING), upper(EPOCH_SELECT_MODE), N_PERM);

%% ====================== MAIN LOOP =======================================
for a = 1:numel(ANALYSIS_STATES)
    analysis_state = ANALYSIS_STATES{a};
    fprintf('\n===== STATE: %s =====\n', analysis_state);

    Data = struct();
    GroupPlot = struct();

    %% ---------- Load and preprocess subject-level data ----------
    for g = 1:numel(GROUPS)
        G = GROUPS{g};
        listing = get_subject_listing(ROOT, METRIC_SUBDIR, G);

        if isempty(listing)
            fprintf('[%s] No files found.\n', G);
            continue;
        end

        fprintf('[%s] %d files found\n', G, numel(listing));

        E_all = [];
        L_all = [];
        SID   = strings(0,1);

        for k = 1:numel(listing)
            subj    = listing(k).subj;
            matPath = listing(k).fullpath;

            [M149, epochStages] = load_sampen_and_scoring( ...
                matPath, SCORING_DIR, G, subj, METRIC_NAME);

            M128 = M149(keep_idx_149to128, :);

            if PRINT_COUNTS
                fprintf('%s\n', epoch_count_string_strict(subj, epochStages, NSEG, EARLY_SEG, LATE_SEG));
            end

            [Esel, Lsel, okSel] = select_E_L_indices_strict( ...
                epochStages, NSEG, EARLY_SEG, LATE_SEG, analysis_state);

            if ~okSel
                fprintf('   [%s] %s excluded: insufficient segmentation.\n', G, subj);
                continue;
            end

            [Esel, Lsel, okClamp] = apply_epoch_clamp( ...
                Esel, Lsel, EPOCH_BALANCING, EPOCH_SELECT_MODE, ...
                MIN_EPOCHS_PER_SEG, RNG_SEED_BALANCING);

            if ~okClamp
                fprintf('   [%s] %s excluded: clamp/min epoch criterion failed.\n', G, subj);
                continue;
            end

            % Aggregate Early/Late epochs per channel
            switch lower(AGG_METHOD)
                case 'mean'
                    vE_raw = mean(M128(:, Esel), 2, 'omitnan')';
                    vL_raw = mean(M128(:, Lsel), 2, 'omitnan')';
                case 'median'
                    vE_raw = median(M128(:, Esel), 2, 'omitnan')';
                    vL_raw = median(M128(:, Lsel), 2, 'omitnan')';
                otherwise
                    error('AGG_METHOD must be ''mean'' or ''median''.');
            end

            % Optional log transform
            switch lower(LOG_STRATEGY)
                case 'post'
                    vE = log(vE_raw + EPS_LOG);
                    vL = log(vL_raw + EPS_LOG);

                case 'pre'
                    logE = log(M128(:, Esel) + EPS_LOG);
                    logL = log(M128(:, Lsel) + EPS_LOG);

                    if strcmpi(AGG_METHOD, 'mean')
                        vE = mean(logE, 2, 'omitnan')';
                        vL = mean(logL, 2, 'omitnan')';
                    else
                        vE = median(logE, 2, 'omitnan')';
                        vL = median(logL, 2, 'omitnan')';
                    end

                case 'none'
                    vE = vE_raw;
                    vL = vL_raw;

                otherwise
                    error('LOG_STRATEGY must be ''post'', ''pre'', or ''none''.');
            end

            E_all(end+1, :) = vE; %#ok<AGROW>
            L_all(end+1, :) = vL; %#ok<AGROW>
            SID(end+1, 1)   = string(subj);
        end

        if isempty(E_all)
            fprintf('   [%s] No valid subjects in %s.\n', G, analysis_state);
            continue;
        end

        Data.(G).E    = E_all;
        Data.(G).L    = L_all;
        Data.(G).Subj = SID;

        GroupPlot.(G).EAR   = mean(E_all, 1, 'omitnan');
        GroupPlot.(G).LAT   = mean(L_all, 1, 'omitnan');
        GroupPlot.(G).nSubj = size(E_all, 1);
    end

    %% ---------- Paired statistics + TFCE ----------
    for g = 1:numel(GROUPS)
        G = GROUPS{g};
        if ~isfield(Data, G), continue; end

        [t_obs, t_perms] = paired_t_per_electrode( ...
            Data.(G).E, Data.(G).L, N_PERM, SIGNFLIP_MODE);

        R = lmeEEG_TFCE_mod(t_obs(:), t_perms, chanlocs, TFCE_EH);
        ES = paired_effect_sizes(Data.(G).E, Data.(G).L);

        Data.(G).within = struct( ...
            'T_Obs',         t_obs(:), ...
            'TFCE_Obs',      R.TFCE_Obs(:), ...
            'P_Values',      R.P_Values(:), ...
            'Mask',          R.Mask(:), ...
            'k_sig',         sum(R.Mask(:)), ...
            'n_per_channel', ES.n(:), ...
            'mean_diff',     ES.mean_diff(:), ...
            'sd_diff',       ES.sd_diff(:), ...
            'cohen_dz',      ES.dz(:), ...
            'hedges_gz',     ES.gz(:));

        fprintf('   [WITHIN TFCE | %s | %s] n=%d, B=%d -> %d/128 significant\n', ...
            analysis_state, G, size(Data.(G).E,1), N_PERM, Data.(G).within.k_sig);

        % Save CSV
        Tout = table(kept128_labels(:), ...
                     Data.(G).within.n_per_channel, ...
                     Data.(G).within.mean_diff, ...
                     Data.(G).within.sd_diff, ...
                     Data.(G).within.cohen_dz, ...
                     Data.(G).within.hedges_gz, ...
                     Data.(G).within.T_Obs, ...
                     Data.(G).within.TFCE_Obs, ...
                     Data.(G).within.P_Values, ...
                     Data.(G).within.Mask, ...
            'VariableNames', {'Electrode','N','MeanDiff','SDDiff', ...
                              'Cohen_dz','Hedges_gz','t_obs', ...
                              'TFCE','p_TFCE','Significant'});

        outCSV = fullfile(OUT_DIR, sprintf( ...
            'WITHIN_TFCE_SampEn_%s_%s_LOG=%s_AGG=%s_CLAMP=%s_%s.csv', ...
            analysis_state, G, LOG_STRATEGY, AGG_METHOD, ...
            EPOCH_BALANCING, EPOCH_SELECT_MODE));
        writetable(Tout, outCSV);
    end

    %% ---------- Figures ----------
    plot_scatter_electrodes_per_group(GroupPlot, GROUPS, analysis_state, OUT_DIR, DPI);
    plot_clean_topogrid_ELD_allgroups(Data, GROUPS, chanlocs, OUT_DIR, analysis_state, DPI);
end

fprintf('\n✅ Done. Output folder: %s\n', OUT_DIR);

%% ========================================================================
%%                                HELPERS
%% ========================================================================

function ensure_dir(d)
if ~exist(d, 'dir')
    mkdir(d);
end
end

function listing = get_subject_listing(ROOT, metricSubdir, G)
p = fullfile(ROOT, metricSubdir, G);
if ~exist(p, 'dir')
    listing = [];
    return;
end

d = dir(fullfile(p, '*_eeg_minimal_*.mat'));
listing = struct('subj', {}, 'fullpath', {}, 'name', {});

for i = 1:numel(d)
    subj = extract_subj_id(d(i).name);
    listing(end+1).subj     = subj; %#ok<AGROW>
    listing(end).fullpath   = fullfile(d(i).folder, d(i).name);
    listing(end).name       = d(i).name;
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

function [M2D, epochStages] = load_sampen_and_scoring(matPath, SCORING_DIR, G, subj, metricVarName)
S = load(matPath);
M2D = pick_sampen_matrix(S, metricVarName);

T = load(fullfile(SCORING_DIR, sprintf('%s_%s_scoring.mat', G, subj)));
epochStages = T.epochStages(:)';

if size(M2D, 2) ~= numel(epochStages)
    error('Epoch mismatch for %s %s.', G, subj);
end
end

function M2D = pick_sampen_matrix(S, preferred)
% Load Sample Entropy matrix from a .mat structure.
% Accepted names: preferred variable, SampleEntropy, SampEn.
% Fallback: first numeric matrix found.

if ~isempty(preferred) && isfield(S, preferred) && isnumeric(S.(preferred))
    M2D = S.(preferred);
    return;
end

candidates = {'SampleEntropy', 'SampEn'};
for i = 1:numel(candidates)
    if isfield(S, candidates{i}) && isnumeric(S.(candidates{i}))
        M2D = S.(candidates{i});
        return;
    end
end

fns = fieldnames(S);
for i = 1:numel(fns)
    V = S.(fns{i});
    if isnumeric(V) && ismatrix(V)
        M2D = V;
        return;
    end
end

error('Sample entropy matrix not found.');
end

function [Esel, Lsel, ok] = select_E_L_indices_strict(epochStages, NSEG, earlySeg, lateSeg, state_sel)
% Build Early/Late epoch indices from NREM segmentation, then keep only
% the requested state.

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
        maskState = (epochStages==2 | epochStages==3);
    case 'N3'
        maskState = (epochStages==3);
    otherwise
        error('state_sel must be ''N2N3'' or ''N3''.');
end

Esel = Eidx(maskState(Eidx));
Lsel = Lidx(maskState(Lidx));

if isempty(Esel) || isempty(Lsel)
    return;
end
ok = true;
end

function [Esel, Lsel, ok] = apply_epoch_clamp(Esel, Lsel, mode, pick, min_each, seed)
% Force Early and Late to the same number of epochs within each subject.

if nargin < 6
    seed = 1234;
end

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
start = floor((N - m)/2) + 1;
ids = pool(start:start+m-1);
end

function [t_obs, t_perms] = paired_t_per_electrode(EARLY, LATE, B, SIGNFLIP_MODE)
% Paired t-statistics per channel with sign-flip permutations.

if nargin < 3 || isempty(B)
    B = 1000;
end
if nargin < 4
    SIGNFLIP_MODE = 'per_channel';
end

[S, C] = size(EARLY);
t_obs   = nan(C,1);
t_perms = nan(B,C);

% Observed statistics
for ch = 1:C
    e = EARLY(:,ch);
    l = LATE(:,ch);
    ok = ~(isnan(e) | isnan(l));
    n = nnz(ok);

    if n < 3
        continue;
    end

    d = l(ok) - e(ok);
    m = mean(d);
    s = std(d, 1);
    t_obs(ch) = m / (s / sqrt(n));
end

% Null distribution
switch lower(SIGNFLIP_MODE)
    case 'per_channel'
        rng(123);
        for b = 1:B
            for ch = 1:C
                e = EARLY(:,ch);
                l = LATE(:,ch);
                ok = ~(isnan(e) | isnan(l));
                n = nnz(ok);

                if n < 3
                    t_perms(b,ch) = NaN;
                    continue;
                end

                sgn = (rand(n,1) > 0.5) * 2 - 1;
                d = (l(ok) - e(ok)) .* sgn;
                m = mean(d);
                s = std(d, 1);
                t_perms(b,ch) = m / (s / sqrt(n));
            end
        end

    case 'per_subject_coherent'
        rng(123);
        sgn_all = ((rand(S,B) > 0.5) * 2 - 1);

        for b = 1:B
            sgn = sgn_all(:,b);
            for ch = 1:C
                e = EARLY(:,ch);
                l = LATE(:,ch);
                ok = ~(isnan(e) | isnan(l));
                n = nnz(ok);

                if n < 3
                    t_perms(b,ch) = NaN;
                    continue;
                end

                d = (l(ok) - e(ok)) .* sgn(ok);
                m = mean(d);
                s = std(d, 1);
                t_perms(b,ch) = m / (s / sqrt(n));
            end
        end

    otherwise
        error('Invalid SIGNFLIP_MODE.');
end
end

function ES = paired_effect_sizes(EARLY, LATE)
% Paired effect sizes per channel:
% dz = mean(diff)/std(diff)
% gz = Hedges-corrected dz

[~, C] = size(EARLY);

n      = nan(1,C);
mDiff  = nan(1,C);
sdDiff = nan(1,C);
dz     = nan(1,C);
gz     = nan(1,C);

for ch = 1:C
    e = EARLY(:,ch);
    l = LATE(:,ch);
    ok = ~(isnan(e) | isnan(l));
    ni = nnz(ok);

    if ni < 3
        continue;
    end

    d = l(ok) - e(ok);
    md = mean(d);
    sd = std(d, 1);

    if sd <= 0
        continue;
    end

    dz_ch = md / sd;
    df = ni - 1;
    J = 1 - 3/(4*df - 1);
    gz_ch = J * dz_ch;

    n(ch)      = ni;
    mDiff(ch)  = md;
    sdDiff(ch) = sd;
    dz(ch)     = dz_ch;
    gz(ch)     = gz_ch;
end

ES = struct('n', n, 'mean_diff', mDiff, 'sd_diff', sdDiff, 'dz', dz, 'gz', gz);
end

function out = epoch_count_string_strict(subj, epochStages, NSEG, earlySeg, lateSeg)
% Print how many NREM and N3 epochs fall in Early/Late windows.

Nepoch = numel(epochStages);
idxNREM = find(epochStages==1 | epochStages==2 | epochStages==3);
nNREM   = numel(idxNREM);

nEarlyNREM = 0;
nLateNREM  = 0;
nEarlyN3   = 0;
nLateN3    = 0;

if numel(idxNREM) >= NSEG
    edges = round(linspace(1, numel(idxNREM)+1, NSEG+1));

    Eidx = [];
    Lidx = [];
    for s = earlySeg
        Eidx = [Eidx, idxNREM(edges(s):edges(s+1)-1)]; %#ok<AGROW>
    end
    for s = lateSeg
        Lidx = [Lidx, idxNREM(edges(s):edges(s+1)-1)]; %#ok<AGROW>
    end

    nEarlyNREM = numel(Eidx);
    nLateNREM  = numel(Lidx);
    nEarlyN3   = sum(epochStages(Eidx)==3);
    nLateN3    = sum(epochStages(Lidx)==3);
end

out = sprintf('  %s: Nepoch=%d | NREM=%d | EARLY(NREM)=%d | LATE(NREM)=%d | EARLY(N3)=%d | LATE(N3)=%d', ...
    subj, Nepoch, nNREM, nEarlyNREM, nLateNREM, nEarlyN3, nLateN3);
end

%% =============================== PLOTTING ================================

function plot_scatter_electrodes_per_group(GroupPlot, GROUPS, analysis_state, OUT_DIR, dpi)
% One figure with one row per group.
% Each panel shows mean Early and Late values across electrodes.

colEarly = [33 113 181]/255;
colLate  = [203 24 29]/255;

allY = [];
for g = 1:numel(GROUPS)
    G = GROUPS{g};
    if ~isfield(GroupPlot, G) || isempty(GroupPlot.(G).EAR), continue; end
    allY = [allY GroupPlot.(G).EAR(:)' GroupPlot.(G).LAT(:)']; %#ok<AGROW>
end

if isempty(allY)
    return;
end

lo  = min(allY,[],'omitnan');
hi  = max(allY,[],'omitnan');
pad = 0.05 * max(hi-lo, eps);

figH = 260 * max(1, numel(GROUPS));
fig  = figure('Color','w','Position',[120 100 1300 figH]);
t    = tiledlayout(numel(GROUPS), 1, 'Padding','compact', 'TileSpacing','compact');

for g = 1:numel(GROUPS)
    G = GROUPS{g};
    ax = nexttile(t); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    if ~isfield(GroupPlot, G) || isempty(GroupPlot.(G).EAR)
        axis off;
        continue;
    end

    e = GroupPlot.(G).EAR(:)';
    l = GroupPlot.(G).LAT(:)';
    x = 1:numel(e);

    scatter(ax, x, e, 22, 'filled', 'MarkerFaceColor', colEarly, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'none');
    scatter(ax, x, l, 22, 'filled', 'MarkerFaceColor', colLate,  'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'none');

    xlim(ax, [1 numel(e)]);
    ylim(ax, [lo-pad hi+pad]);
    ylabel(ax, G, 'FontWeight', 'bold');

    if g == numel(GROUPS)
        xlabel(ax, 'Electrode (1..128)');
    else
        set(ax, 'XTickLabel', []);
    end
end

title(t, sprintf('Sample Entropy — %s — Early (blue) vs Late (red)', analysis_state), 'FontWeight', 'bold');

outPNG = fullfile(OUT_DIR, sprintf('ALL_SampEn_%s_SCATTER_ELECTRODES_PERGROUP.png', analysis_state));
exportgraphics(fig, outPNG, 'Resolution', dpi);
fprintf('Saved: %s\n', outPNG);
end

function plot_clean_topogrid_ELD_allgroups(Data, GROUPS, chanlocs, OUT_DIR, analysis_state, dpi)
% One figure:
%   rows = groups
%   cols = Early / Late / Diff
% Diff maps include TFCE-significant electrodes.

nG   = numel(GROUPS);
have = false(1,nG);
Emap = [];
Lmap = [];
Dmap = [];
SigMask = cell(1,nG);

for g = 1:nG
    G = GROUPS{g};
    if isfield(Data, G) && isfield(Data.(G), 'E') && ~isempty(Data.(G).E)
        eM = mean(Data.(G).E, 1, 'omitnan');
        lM = mean(Data.(G).L, 1, 'omitnan');
        dM = lM - eM;

        Emap = [Emap; eM]; %#ok<AGROW>
        Lmap = [Lmap; lM]; %#ok<AGROW>
        Dmap = [Dmap; dM]; %#ok<AGROW>
        have(g) = true;

        if isfield(Data.(G), 'within') && isfield(Data.(G).within, 'Mask') && ~isempty(Data.(G).within.Mask)
            SigMask{g} = logical(Data.(G).within.Mask(:));
        else
            SigMask{g} = false(numel(eM),1);
        end
    end
end

if ~any(have)
    return;
end

allEL  = [Emap; Lmap];
climEL = [min(allEL,[],'all','omitnan'), max(allEL,[],'all','omitnan')];

maxAbsD = max(abs(Dmap),[],'all','omitnan');
climD   = [-maxAbsD maxAbsD];

figH = 260 * nnz(have) + 80;
figW = 3 * 260 + 80;
fig  = figure('Color','w','Position',[100 100 figW figH]);
t    = tiledlayout(nnz(have), 3, 'Padding','compact', 'TileSpacing','compact');

rowIdx = 0;
for g = 1:nG
    G = GROUPS{g};
    if ~have(g), continue; end
    rowIdx = rowIdx + 1;

    eM = mean(Data.(G).E, 1, 'omitnan');
    lM = mean(Data.(G).L, 1, 'omitnan');
    dM = lM - eM;

    % Early
    ax1 = nexttile(t, (rowIdx-1)*3 + 1);
    axes(ax1); cla(ax1); axis(ax1,'off');
    draw_clean_topo(eM(:)', chanlocs, climEL);
    caxis(ax1, climEL);
    colorbar('Location','eastoutside');
    ylabel(ax1, G, 'FontWeight', 'bold');

    % Late
    ax2 = nexttile(t, (rowIdx-1)*3 + 2);
    axes(ax2); cla(ax2); axis(ax2,'off');
    draw_clean_topo(lM(:)', chanlocs, climEL);
    caxis(ax2, climEL);
    colorbar('Location','eastoutside');

    % Diff
    ax3 = nexttile(t, (rowIdx-1)*3 + 3);
    axes(ax3); cla(ax3); axis(ax3,'off');
    draw_diff_topo_clean_with_sig(dM(:)', chanlocs, climD, SigMask{g});
    caxis(ax3, climD);
    colorbar('Location','eastoutside');
end

outPNG = fullfile(OUT_DIR, sprintf('ALL_SampEn_%s_TOPOGRID_CLEAN_ELD_SIG.png', analysis_state));
exportgraphics(fig, outPNG, 'Resolution', dpi);
fprintf('Saved: %s\n', outPNG);
end

function draw_clean_topo(topo, chanlocs, clim)
if exist('topoplot','file') == 2
    try
        topoplot(topo(:), chanlocs, 'maplimits', clim, 'electrodes', 'off', 'style', 'both');
        axis tight; axis off;
        caxis(clim);
        return;
    catch
    end
end

xs = [chanlocs.X];
ys = [chanlocs.Y];
scatter(xs, ys, 80, topo, 'filled');
axis equal; axis tight; axis off;
caxis(clim);
end

function draw_diff_topo_clean_with_sig(topo, chanlocs, clim, SigMask)
if nargin < 4 || isempty(SigMask)
    SigMask = false(numel(topo),1);
else
    SigMask = logical(SigMask(:));
end

if exist('topoplot','file') == 2
    try
        sigCh = find(SigMask);
        if isempty(sigCh)
            topoplot(topo(:), chanlocs, 'maplimits', clim, 'electrodes', 'off', 'style', 'both');
        else
            topoplot(topo(:), chanlocs, 'maplimits', clim, 'electrodes', 'off', 'style', 'both', ...
                'emarker2', {sigCh, 'o', [1 1 1], 3, 1.5});
        end
        axis tight; axis off;
        caxis(clim);
        return;
    catch
    end
end

xs = [chanlocs.X];
ys = [chanlocs.Y];
scatter(xs, ys, 80, topo, 'filled');
axis equal; axis tight; axis off;
caxis(clim);

if any(SigMask)
    hold on;
    idx = find(SigMask);
    plot([chanlocs(idx).X], [chanlocs(idx).Y], 'wo', 'MarkerSize', 6, 'LineWidth', 1.5);
    plot([chanlocs(idx).X], [chanlocs(idx).Y], 'ko', 'MarkerSize', 6, 'LineWidth', 0.8);
end
end