%% ========================================================================
% SWA — EARLY vs LATE (WITHIN) — TFCE FROM SWA_CLAMP OUTPUT
%
% DESCRIPTION:
%   This script performs within-group Early vs Late comparisons using
%   SWA_CLAMP output that was already matched/clamped at the subject level.
%
%   Only two metrics are analyzed:
%       - classic
%       - density
%
%   For each group and metric:
%       1) Early and Late matrices are loaded from globalRes
%       2) Optional log-transform is applied at the subject level
%       3) Paired t-statistics are computed electrode-wise
%       4) Sign-flip permutations are used to build the null distribution
%       5) TFCE is applied using lmeEEG_TFCE_mod
%       6) CSV summaries and figures are exported
%
% INPUT:
%   - ...\for_reviewers\SWA_CLAMP\globalRes_SWA_CLAMP.mat
%
% OUTPUT:
%   - One CSV per group × metric
%   - One topogrid figure per metric
%   - One overall paired figure per metric
%
% REQUIREMENTS:
%   - EEGLAB (topoplot) for topographic plots
%   - lmeEEG_TFCE_mod on MATLAB path
% ========================================================================

clear; clc; close all; rng(42);

%% ====================== USER SETTINGS ===================================
% Input
ROOT_SWA = '...\for_reviewers\SWA_CLAMP';
IN_MAT   = fullfile(ROOT_SWA, 'globalRes_SWA_CLAMP.mat');

% Output
OUTROOT  = '...\Annals_entropy';
OUT_DIR  = fullfile(OUTROOT, 'GroupAnalysis', 'WITHIN_SWA_FROM_CLAMP');
ensure_dir(OUT_DIR);

% Groups available in globalRes
GROUPS = {'CTL','ADV','DYS'};

% Band to analyze
BAND = [0.5 4];

% Only these two metrics are analyzed
METRICS = {'classic','density'};

% Log-transform on subject-level values
% Both classic and density are positive, so log-transform is allowed
LOG_STRATEGY = 'post';   % 'post' or 'none'
EPS_LOG      = 1e-8;

% TFCE / permutations
TFCE_EH       = [0.66 2];
N_PERM        = 1000;
SIGNFLIP_MODE = 'per_channel';   % 'per_channel' or 'per_subject_coherent'

% Plot options
DPI = 300;
PLOT_TOPOGRID     = true;
PLOT_OVERALL_PAIRED = true;

%% ====================== LOAD DATA =======================================
if ~exist(IN_MAT, 'file')
    error('File not found: %s', IN_MAT);
end

S = load(IN_MAT, 'globalRes');
globalRes = S.globalRes;

% Print metadata if present
if isfield(globalRes, 'meta')
    M = globalRes.meta;
    fprintf('\n=== SWA FROM CLAMP ===\n');
    fprintf('State: %s | NSEG=%d | early=%s | late=%s | epoch_len=%ds\n', ...
        string(M.STATE_SEL), M.NSEG, mat2str(M.earlySeg), mat2str(M.lateSeg), M.epoch_len_s);
    fprintf('Clamp: %s / %s | min=%.1f min (=%d epochs)\n', ...
        string(M.EPOCH_BALANCING), string(M.EPOCH_SELECT_MODE), ...
        M.MIN_MATCH_MINUTES_N3, M.min_epochs);
end

% Build band tag used in globalRes
btag = band_tag(BAND);

% Check required functions
if exist('lmeEEG_TFCE_mod', 'file') ~= 2
    error('lmeEEG_TFCE_mod not found on MATLAB path.');
end

if exist('topoplot', 'file') ~= 2
    warning('topoplot not found. Topographic plots will fall back to scatter mode.');
end

fprintf('\nUsing BAND = [%.1f %.1f] | tag = %s | LOG = %s | N_PERM = %d | SIGNFLIP = %s\n', ...
    BAND(1), BAND(2), btag, upper(LOG_STRATEGY), N_PERM, SIGNFLIP_MODE);

%% ====================== MAIN ANALYSIS LOOP ==============================
for im = 1:numel(METRICS)
    METRIC = METRICS{im};
    fprintf('\n============================================================\n');
    fprintf('METRIC: %s\n', upper(METRIC));
    fprintf('============================================================\n');

    Data = struct();
    kept_labels = [];

    %% ---------- Load Early/Late matrices for all groups ----------
    for g = 1:numel(GROUPS)
        G = GROUPS{g};

        if ~isfield(globalRes, G)
            fprintf('[%s] not found in globalRes -> skip\n', G);
            continue;
        end
        if ~isfield(globalRes.(G), btag)
            fprintf('[%s] band %s not found -> skip\n', G, btag);
            continue;
        end
        if ~isfield(globalRes.(G).(btag).early, METRIC)
            fprintf('[%s] metric %s not found -> skip\n', G, METRIC);
            continue;
        end

        % Early/Late matrices come as [nCh x nSubj]
        E = globalRes.(G).(btag).early.(METRIC);
        L = globalRes.(G).(btag).late.(METRIC);

        % Convert to [nSubj x nCh]
        E = E';
        L = L';

        % Optional log-transform
        switch lower(LOG_STRATEGY)
            case 'post'
                E = log(E + EPS_LOG);
                L = log(L + EPS_LOG);
            case 'none'
                % no operation
            otherwise
                error('LOG_STRATEGY must be ''post'' or ''none''.');
        end

        Data.(G).E = E;
        Data.(G).L = L;
        Data.(G).Subj = globalRes.(G).subj(:);
        Data.(G).chanlocs = globalRes.(G).chanlocs;

        if isempty(kept_labels) && isfield(globalRes.(G), 'labels')
            kept_labels = string(globalRes.(G).labels(:));
        end

        fprintf('[%s] loaded: %d subjects, %d channels\n', G, size(E,1), size(E,2));
    end

    %% ---------- Within-group TFCE analysis ----------
    GroupPlot = struct();

    for g = 1:numel(GROUPS)
        G = GROUPS{g};
        if ~isfield(Data, G), continue; end

        E = Data.(G).E;
        L = Data.(G).L;
        chanlocs = Data.(G).chanlocs;

        % Mean Early/Late maps across subjects
        GroupPlot.(G).EAR  = mean(E, 1, 'omitnan');
        GroupPlot.(G).LAT  = mean(L, 1, 'omitnan');
        GroupPlot.(G).DIFF = mean(L - E, 1, 'omitnan');
        GroupPlot.(G).nSubj = size(E,1);

        % Electrode-wise paired t + permutations
        [t_obs, t_perms] = paired_t_per_electrode(E, L, N_PERM, SIGNFLIP_MODE);

        % TFCE
        Rtfce = lmeEEG_TFCE_mod(t_obs(:), t_perms, chanlocs, TFCE_EH);

        % Paired effect sizes
        ES = paired_effect_sizes(E, L);

        Data.(G).within = struct( ...
            'T_Obs',    t_obs(:), ...
            'TFCE_Obs', Rtfce.TFCE_Obs(:), ...
            'P_Values', Rtfce.P_Values(:), ...
            'Mask',     Rtfce.Mask(:), ...
            'k_sig',    sum(Rtfce.Mask(:)), ...
            'n_per_channel', ES.n(:), ...
            'mean_diff', ES.mean_diff(:), ...
            'sd_diff',   ES.sd_diff(:), ...
            'cohen_dz',  ES.dz(:), ...
            'hedges_gz', ES.gz(:));

        fprintf('[WITHIN TFCE | %s | %s] nSubj=%d, nPerm=%d -> %d/%d significant electrodes\n', ...
            G, upper(METRIC), size(E,1), N_PERM, Data.(G).within.k_sig, size(E,2));

        % Save CSV
        Electrode = kept_labels;
        if isempty(Electrode)
            Electrode = string((1:size(E,2))');
        end

        Tout = table(Electrode(:), ...
            Data.(G).within.n_per_channel, ...
            Data.(G).within.mean_diff, ...
            Data.(G).within.sd_diff, ...
            Data.(G).within.cohen_dz, ...
            Data.(G).within.hedges_gz, ...
            Data.(G).within.T_Obs, ...
            Data.(G).within.TFCE_Obs, ...
            Data.(G).within.P_Values, ...
            Data.(G).within.Mask, ...
            'VariableNames', {'Electrode','N','MeanDiff','SDDiff','Cohen_dz','Hedges_gz','t_obs','TFCE','p_TFCE','Significant'});

        outCSV = fullfile(OUT_DIR, sprintf('WITHIN_TFCE_SWA_%s_%s_BAND=%s_LOG=%s_B=%d.csv', ...
            upper(METRIC), G, btag, upper(LOG_STRATEGY), N_PERM));
        writetable(Tout, outCSV);
    end

    %% ---------- Figures ----------
    TITLE_SFX = sprintf('SWA %s | band=[%.1f %.1f] | LOG=%s', ...
        upper(METRIC), BAND(1), BAND(2), upper(LOG_STRATEGY));

    if PLOT_TOPOGRID
        plot_topogrid_ELD_allgroups(Data, GROUPS, OUT_DIR, TITLE_SFX, DPI);
    end

    if PLOT_OVERALL_PAIRED
        plot_overall_paired(Data, GROUPS, TITLE_SFX, OUT_DIR, DPI);
    end
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

function tag = band_tag(fb)
tag = sprintf('B_%g_%g', fb(1), fb(2));
tag = strrep(tag, '.', 'p');
end

function [t_obs, t_perms] = paired_t_per_electrode(EARLY, LATE, B, SIGNFLIP_MODE)
% Compute paired t-statistics electrode-wise and build a sign-flip null.
%
% INPUT:
%   EARLY, LATE : [nSubj x nCh]
%   B           : number of permutations
%   SIGNFLIP_MODE:
%       'per_channel'          -> independent sign flips for each channel
%       'per_subject_coherent' -> same sign pattern across channels per subject

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

% Null distribution via sign-flip
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
% Compute paired effect sizes per channel:
%   dz = mean(diff) / std(diff)
%   gz = Hedges-corrected dz

[~, C] = size(EARLY);

n  = nan(1,C);
md = nan(1,C);
sd = nan(1,C);
dz = nan(1,C);
gz = nan(1,C);

for ch = 1:C
    e = EARLY(:,ch);
    l = LATE(:,ch);
    ok = ~(isnan(e) | isnan(l));

    ni = nnz(ok);
    if ni < 3
        continue;
    end

    d = l(ok) - e(ok);
    m = mean(d);
    s = std(d, 1);

    if s <= 0
        continue;
    end

    dz_ch = m / s;
    df = ni - 1;

    if df > 0
        J = 1 - 3 / (4*df - 1);
    else
        J = NaN;
    end

    gz_ch = J * dz_ch;

    n(ch)  = ni;
    md(ch) = m;
    sd(ch) = s;
    dz(ch) = dz_ch;
    gz(ch) = gz_ch;
end

ES = struct('n', n, 'mean_diff', md, 'sd_diff', sd, 'dz', dz, 'gz', gz);
end

function plot_topogrid_ELD_allgroups(Data, GROUPS, OUT_DIR, title_sfx, dpi)
% Plot a topographic grid:
%   rows = groups
%   cols = Early / Late / Diff
% Diff maps include TFCE significant electrodes.

Emap = [];
Lmap = [];
Dmap = [];
have = false(1, numel(GROUPS));
SigMask = cell(1, numel(GROUPS));
chanlocs_ref = [];

for g = 1:numel(GROUPS)
    G = GROUPS{g};
    if ~isfield(Data, G), continue; end

    E = Data.(G).E;
    L = Data.(G).L;

    eM = mean(E,1,'omitnan');
    lM = mean(L,1,'omitnan');
    dM = lM - eM;

    Emap = [Emap; eM]; %#ok<AGROW>
    Lmap = [Lmap; lM]; %#ok<AGROW>
    Dmap = [Dmap; dM]; %#ok<AGROW>
    have(g) = true;

    if isfield(Data.(G), 'within') && isfield(Data.(G).within, 'Mask')
        SigMask{g} = logical(Data.(G).within.Mask(:));
    else
        SigMask{g} = false(numel(eM),1);
    end

    chanlocs_ref = Data.(G).chanlocs;
end

if ~any(have)
    return;
end

% Common color limits
allEL  = [Emap; Lmap];
climEL = [min(allEL,[],'all','omitnan'), max(allEL,[],'all','omitnan')];

maxAbsD = max(abs(Dmap),[],'all','omitnan');
climD   = [-maxAbsD maxAbsD];

nRows = nnz(have);
fig = figure('Color','w','Position',[100 100 980 260*nRows+60]);
t = tiledlayout(nRows, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

rowIdx = 0;
for g = 1:numel(GROUPS)
    G = GROUPS{g};
    if ~have(g), continue; end
    rowIdx = rowIdx + 1;

    E = Data.(G).E;
    L = Data.(G).L;

    eM = mean(E,1,'omitnan');
    lM = mean(L,1,'omitnan');
    dM = lM - eM;

    % Early
    ax1 = nexttile(t, (rowIdx-1)*3 + 1);
    axes(ax1); cla(ax1); axis(ax1,'off');
    draw_topo(eM(:), chanlocs_ref, climEL);
    caxis(ax1, climEL);
    cb1 = colorbar('Location','eastoutside');
    cb1.FontWeight = 'bold';
    ylabel(ax1, G, 'FontWeight', 'bold');

    % Late
    ax2 = nexttile(t, (rowIdx-1)*3 + 2);
    axes(ax2); cla(ax2); axis(ax2,'off');
    draw_topo(lM(:), chanlocs_ref, climEL);
    caxis(ax2, climEL);
    cb2 = colorbar('Location','eastoutside');
    cb2.FontWeight = 'bold';

    % Diff
    ax3 = nexttile(t, (rowIdx-1)*3 + 3);
    axes(ax3); cla(ax3); axis(ax3,'off');
    draw_topo_diff_with_sig(dM(:), chanlocs_ref, climD, SigMask{g});
    caxis(ax3, climD);
    cb3 = colorbar('Location','eastoutside');
    cb3.FontWeight = 'bold';
end

outPNG = fullfile(OUT_DIR, sprintf('SWA_%s_TOPOGRID_ELD_SIG.png', regexprep(title_sfx,'[^A-Za-z0-9]','_')));
exportgraphics(fig, outPNG, 'Resolution', dpi);
fprintf('Saved: %s\n', outPNG);
end

function draw_topo(vals, chanlocs, clim)
if exist('topoplot', 'file') == 2
    topoplot(vals, chanlocs, 'electrodes', 'off', 'style', 'both', 'maplimits', clim);
    axis tight; axis off; caxis(clim);
else
    xs = [chanlocs.X];
    ys = [chanlocs.Y];
    scatter(xs, ys, 80, vals, 'filled');
    axis equal; axis off; caxis(clim);
end
end

function draw_topo_diff_with_sig(vals, chanlocs, clim, SigMask)
SigMask = logical(SigMask(:));

if exist('topoplot', 'file') == 2
    sigCh = find(SigMask);
    if isempty(sigCh)
        topoplot(vals, chanlocs, 'electrodes', 'off', 'style', 'both', 'maplimits', clim);
    else
        topoplot(vals, chanlocs, 'electrodes', 'off', 'style', 'both', 'maplimits', clim, ...
            'emarker2', {sigCh, 'o', [1 1 1], 3, 1.5});
    end
    axis tight; axis off; caxis(clim);
else
    xs = [chanlocs.X];
    ys = [chanlocs.Y];
    scatter(xs, ys, 80, vals, 'filled');
    axis equal; axis off; caxis(clim);
    if any(SigMask)
        hold on;
        idx = find(SigMask);
        plot([chanlocs(idx).X], [chanlocs(idx).Y], 'wo', 'MarkerSize', 6, 'LineWidth', 1.5);
        plot([chanlocs(idx).X], [chanlocs(idx).Y], 'ko', 'MarkerSize', 6, 'LineWidth', 0.8);
    end
end
end

function plot_overall_paired(Data, GROUPS, title_sfx, OUT_DIR, dpi)
% Plot subject-level Early vs Late values after averaging across electrodes.

fig = figure('Color','w','Position',[100 100 1200 420]);
t = tiledlayout(1, numel(GROUPS), 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(sprintf('%s | Overall paired Early vs Late', title_sfx), 'FontWeight', 'bold');

for g = 1:numel(GROUPS)
    G = GROUPS{g};
    ax = nexttile(t); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    if ~isfield(Data, G)
        axis off;
        continue;
    end

    e = mean(Data.(G).E, 2, 'omitnan');
    l = mean(Data.(G).L, 2, 'omitnan');
    ok = isfinite(e) & isfinite(l);
    e = e(ok);
    l = l(ok);

    for i = 1:numel(e)
        plot([1 2], [e(i) l(i)], '-o', 'LineWidth', 1);
    end

    xlim([0.8 2.2]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Early','Late'});

    if numel(e) >= 3
        p = signrank(e, l);
        title(sprintf('%s | p=%s | n=%d', G, fmt_p(p), numel(e)));
    else
        title(sprintf('%s | n=%d', G, numel(e)));
    end
end

outPNG = fullfile(OUT_DIR, sprintf('SWA_%s_OVERALL_PAIRED.png', regexprep(title_sfx,'[^A-Za-z0-9]','_')));
exportgraphics(fig, outPNG, 'Resolution', dpi);
fprintf('Saved: %s\n', outPNG);
end

function s = fmt_p(p)
if ~isfinite(p)
    s = 'NA';
elseif p < 1e-4
    s = sprintf('%.1e', p);
else
    s = sprintf('%.4f', p);
end
end