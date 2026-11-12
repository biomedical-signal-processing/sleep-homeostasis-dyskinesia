%% ========================================================================
% BASIC QC SUMMARY FROM EPOCHED .SET FILES
%
% DESCRIPTION:
% This script extracts basic QC information directly from epoched EEGLAB
% .set files, assuming that the relevant epoch-level information is already
% stored inside each dataset.
%
% The script does NOT reconstruct stages/arousals from sample-wise signals
% and does NOT match raw vs processed datasets. It simply reads each .set
% file as-is and summarizes the available QC fields.
%
% EXPECTED OPTIONAL FIELDS IN EACH .set:
%   - EEG.swa_scoring.epochStages
%       or EEG.etc.swa_scoring.epochStages
%       -> epoch-wise sleep stage labels
%
%   - EEG.swa_scoring.epochArousals
%       or EEG.etc.swa_scoring.epochArousals
%       -> epoch-wise arousal flags (0/1)
%
%   - EEG.csc_marked_trials
%       or EEG.etc.csc_marked_trials
%       -> epoch-wise visually rejected epochs (0/1)
%
%   - EEG.chaninfo.icachansind
%       -> channels used for ICA
%
%   - EEG.etc.ic_classification.ICLabel.classifications
%       -> ICLabel probabilities for automatic bad IC count
%
% OUTPUTS:
%   - QC_basic_PerSubject.csv
%   - QC_basic_GroupSummary.csv
%   - FIG_ArousalPct_GroupBar.png         (if arousal info exists)
%   - FIG_VisualBadPct_GroupBar.png       (if visual rejection exists)
%   - TABLE_QC_BasicRecap.csv
%
% REQUIREMENTS:
%   - EEGLAB (pop_loadset)
% ========================================================================

clear; clc; close all;

%% ===================== USER SETTINGS =====================
dataRoot = '...\Annals_processed';
groups   = {'CTL','ADV','DYS'};

% Sleep stage codes
ST_WAKE = 0;
ST_N1   = 1;
ST_N2   = 2;
ST_N3   = 3;
ST_REM  = 5;

% Expected scalp channels in your pipeline
N_SCALP = 149;

outDir = fullfile(dataRoot, '_QC_BASIC');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ===================== MAIN LOOP =====================
rows = struct([]);

for ig = 1:numel(groups)
    G = groups{ig};

    files = dir(fullfile(dataRoot, G, '**', '*.set'));
    fprintf('\n[%s] .set files found: %d\n', G, numel(files));

    for k = 1:numel(files)
        fpath = fullfile(files(k).folder, files(k).name);
        sid   = erase(files(k).name, '.set');

        fprintf('  -> %s / %s\n', G, sid);

        try
            EEG = pop_loadset('filename', files(k).name, 'filepath', files(k).folder);
        catch ME
            warning('Could not load %s: %s', files(k).name, ME.message);
            continue;
        end

        %% ---------- Basic dataset info ----------
        sr      = double(EEG.srate);
        nEpoch  = EEG.trials;
        nChan   = EEG.nbchan;

        %% ---------- Read epoch-wise stage information ----------
        epochStages = [];
        if isfield(EEG, 'swa_scoring') && isfield(EEG.swa_scoring, 'epochStages')
            epochStages = EEG.swa_scoring.epochStages(:);
        elseif isfield(EEG, 'etc') && isfield(EEG.etc, 'swa_scoring') && isfield(EEG.etc.swa_scoring, 'epochStages')
            epochStages = EEG.etc.swa_scoring.epochStages(:);
        end

        %% ---------- Read epoch-wise arousal information ----------
        epochArousals = [];
        if isfield(EEG, 'swa_scoring') && isfield(EEG.swa_scoring, 'epochArousals')
            epochArousals = EEG.swa_scoring.epochArousals(:);
        elseif isfield(EEG, 'etc') && isfield(EEG.etc, 'swa_scoring') && isfield(EEG.etc.swa_scoring, 'epochArousals')
            epochArousals = EEG.etc.swa_scoring.epochArousals(:);
        end

        %% ---------- Read visual rejection information ----------
        csc = [];
        if isfield(EEG, 'etc') && isfield(EEG.etc, 'csc_marked_trials') && ~isempty(EEG.etc.csc_marked_trials)
            csc = EEG.etc.csc_marked_trials(:);
        elseif isfield(EEG, 'csc_marked_trials') && ~isempty(EEG.csc_marked_trials)
            csc = EEG.csc_marked_trials(:);
        end

        %% ---------- Align lengths if needed ----------
        nUse = nEpoch;

        if ~isempty(epochStages)
            nUse = min(nUse, numel(epochStages));
        end
        if ~isempty(epochArousals)
            nUse = min(nUse, numel(epochArousals));
        end
        if ~isempty(csc)
            nUse = min(nUse, numel(csc));
        end

        if ~isempty(epochStages),   epochStages   = epochStages(1:nUse); end
        if ~isempty(epochArousals), epochArousals = logical(epochArousals(1:nUse)); end
        if ~isempty(csc),           csc           = logical(csc(1:nUse)); end

        %% ---------- Stage counts ----------
        nEpoch_WAKE = NaN;
        nEpoch_REM  = NaN;
        nEpoch_N1   = NaN;
        nEpoch_N2   = NaN;
        nEpoch_N3   = NaN;
        nEpoch_NREM = NaN;

        isWake = [];
        isREM  = [];
        isN1   = [];
        isN2   = [];
        isN3   = [];
        isNREM = [];

        if ~isempty(epochStages)
            isWake = (epochStages == ST_WAKE);
            isREM  = (epochStages == ST_REM);
            isN1   = (epochStages == ST_N1);
            isN2   = (epochStages == ST_N2);
            isN3   = (epochStages == ST_N3);
            isNREM = (isN1 | isN2 | isN3);

            nEpoch_WAKE = sum(isWake);
            nEpoch_REM  = sum(isREM);
            nEpoch_N1   = sum(isN1);
            nEpoch_N2   = sum(isN2);
            nEpoch_N3   = sum(isN3);
            nEpoch_NREM = sum(isNREM);
        end

        %% ---------- Arousal-related QC ----------
        ArousalBad_total = NaN;
        ArousalBad_NREM  = NaN;
        ArousalBad_N3    = NaN;

        ArousalBadPct_total = NaN;
        ArousalBadPct_NREM  = NaN;
        ArousalBadPct_N3    = NaN;

        if ~isempty(epochArousals)
            ArousalBad_total = sum(epochArousals);
            ArousalBadPct_total = ArousalBad_total / max(nUse, 1);

            if ~isempty(epochStages)
                ArousalBad_NREM = sum(epochArousals & isNREM);
                ArousalBad_N3   = sum(epochArousals & isN3);

                ArousalBadPct_NREM = ArousalBad_NREM / max(nEpoch_NREM, 1);
                ArousalBadPct_N3   = ArousalBad_N3   / max(nEpoch_N3, 1);
            end
        end

        %% ---------- Visual rejection QC ----------
        VisualBad_total = NaN;
        VisualBad_WAKE  = NaN;
        VisualBad_REM   = NaN;
        VisualBad_N1    = NaN;
        VisualBad_N2    = NaN;
        VisualBad_N3    = NaN;
        VisualBad_NREM  = NaN;

        VisualBadPct_total = NaN;
        VisualBadPct_WAKE  = NaN;
        VisualBadPct_REM   = NaN;
        VisualBadPct_N1    = NaN;
        VisualBadPct_N2    = NaN;
        VisualBadPct_N3    = NaN;
        VisualBadPct_NREM  = NaN;

        if ~isempty(csc)
            VisualBad_total = sum(csc);
            VisualBadPct_total = VisualBad_total / max(nUse, 1);

            if ~isempty(epochStages)
                VisualBad_WAKE = sum(csc & isWake);
                VisualBad_REM  = sum(csc & isREM);
                VisualBad_N1   = sum(csc & isN1);
                VisualBad_N2   = sum(csc & isN2);
                VisualBad_N3   = sum(csc & isN3);
                VisualBad_NREM = sum(csc & isNREM);

                VisualBadPct_WAKE = VisualBad_WAKE / max(nEpoch_WAKE, 1);
                VisualBadPct_REM  = VisualBad_REM  / max(nEpoch_REM, 1);
                VisualBadPct_N1   = VisualBad_N1   / max(nEpoch_N1, 1);
                VisualBadPct_N2   = VisualBad_N2   / max(nEpoch_N2, 1);
                VisualBadPct_N3   = VisualBad_N3   / max(nEpoch_N3, 1);
                VisualBadPct_NREM = VisualBad_NREM / max(nEpoch_NREM, 1);
            end
        end

        %% ---------- ICA QC ----------
        ICA_nCh_used         = NaN;
        ICA_nCh_interpolated = NaN;
        ICA_nIC_bad_auto     = NaN;

        if isfield(EEG, 'chaninfo') && isfield(EEG.chaninfo, 'icachansind') && ~isempty(EEG.chaninfo.icachansind)
            ICA_nCh_used = numel(EEG.chaninfo.icachansind);
            ICA_nCh_interpolated = N_SCALP - ICA_nCh_used;
        end

        if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') && ...
           isfield(EEG.etc.ic_classification, 'ICLabel') && ...
           isfield(EEG.etc.ic_classification.ICLabel, 'classifications') && ...
           ~isempty(EEG.etc.ic_classification.ICLabel.classifications)

            C = EEG.etc.ic_classification.ICLabel.classifications;

            badICs_auto = unique([
                find(C(:,3) > 0.8);  % eyes
                find(C(:,2) > 0.8);  % muscle
                find(C(:,4) > 0.8);  % heart
                find(C(:,6) > 0.7);  % channel noise
                find(C(:,7) > 0.9)   % other
            ]);

            ICA_nIC_bad_auto = numel(badICs_auto);
        end

        %% ---------- Save one row ----------
        r = numel(rows) + 1;

        rows(r).Group   = string(G);
        rows(r).Subject = string(sid);
        rows(r).FilePath = string(fpath);

        rows(r).srate   = sr;
        rows(r).nChan   = nChan;
        rows(r).nEpoch_total = nUse;

        rows(r).nEpoch_WAKE = nEpoch_WAKE;
        rows(r).nEpoch_REM  = nEpoch_REM;
        rows(r).nEpoch_N1   = nEpoch_N1;
        rows(r).nEpoch_N2   = nEpoch_N2;
        rows(r).nEpoch_N3   = nEpoch_N3;
        rows(r).nEpoch_NREM = nEpoch_NREM;

        rows(r).ArousalBad_total    = ArousalBad_total;
        rows(r).ArousalBad_NREM     = ArousalBad_NREM;
        rows(r).ArousalBad_N3       = ArousalBad_N3;
        rows(r).ArousalBadPct_total = ArousalBadPct_total;
        rows(r).ArousalBadPct_NREM  = ArousalBadPct_NREM;
        rows(r).ArousalBadPct_N3    = ArousalBadPct_N3;

        rows(r).VisualBad_total     = VisualBad_total;
        rows(r).VisualBad_WAKE      = VisualBad_WAKE;
        rows(r).VisualBad_REM       = VisualBad_REM;
        rows(r).VisualBad_N1        = VisualBad_N1;
        rows(r).VisualBad_N2        = VisualBad_N2;
        rows(r).VisualBad_N3        = VisualBad_N3;
        rows(r).VisualBad_NREM      = VisualBad_NREM;
        rows(r).VisualBadPct_total  = VisualBadPct_total;
        rows(r).VisualBadPct_WAKE   = VisualBadPct_WAKE;
        rows(r).VisualBadPct_REM    = VisualBadPct_REM;
        rows(r).VisualBadPct_N1     = VisualBadPct_N1;
        rows(r).VisualBadPct_N2     = VisualBadPct_N2;
        rows(r).VisualBadPct_N3     = VisualBadPct_N3;
        rows(r).VisualBadPct_NREM   = VisualBadPct_NREM;

        rows(r).ICA_nCh_used         = ICA_nCh_used;
        rows(r).ICA_nCh_interpolated = ICA_nCh_interpolated;
        rows(r).ICA_nIC_bad_auto     = ICA_nIC_bad_auto;
    end
end

%% ===================== BUILD TABLES =====================
Tsubj = struct2table(rows);

csv1 = fullfile(outDir, 'QC_basic_PerSubject.csv');
writetable(Tsubj, csv1);

% Group summary: mean ± SD for numeric variables
vars = Tsubj.Properties.VariableNames;
isNum = false(size(vars));
for i = 1:numel(vars)
    isNum(i) = isnumeric(Tsubj.(vars{i}));
end
numVars = vars(isNum);

out = struct([]);
for ig = 1:numel(groups)
    G = string(groups{ig});
    idx = (categorical(Tsubj.Group, groups, groups) == G);

    S = Tsubj(idx, :);
    rr = struct();
    rr.Group = G;
    rr.N = height(S);

    for v = 1:numel(numVars)
        vn = numVars{v};
        x = S.(vn);
        rr.([vn '_mean']) = mean(x, 'omitnan');
        rr.([vn '_sd'])   = std(x, 0, 'omitnan');
    end

    out = [out; rr]; %#ok<AGROW>
end

Tgrp = struct2table(out);
csv2 = fullfile(outDir, 'QC_basic_GroupSummary.csv');
writetable(Tgrp, csv2);

fprintf('\nDONE.\nSaved:\n  %s\n  %s\n', csv1, csv2);

%% ===================== SIMPLE PLOTS =====================
grp = categorical(Tsubj.Group, groups, groups);

% ---------- Figure 1: Arousal percentages ----------
if any(~isnan(Tsubj.ArousalBadPct_NREM)) || any(~isnan(Tsubj.ArousalBadPct_N3))

    metricsA = {'ArousalBadPct_NREM','ArousalBadPct_N3'};
    labelsA  = {'NREM','N3'};

    M = nan(numel(groups), numel(metricsA));
    S = nan(numel(groups), numel(metricsA));

    for ig = 1:numel(groups)
        idx = (grp == groups{ig});
        for im = 1:numel(metricsA)
            x = Tsubj.(metricsA{im})(idx);
            M(ig,im) = mean(x, 'omitnan');
            S(ig,im) = std(x, 0, 'omitnan');
        end
    end

    figure('Color','w');
    b = bar(M, 'grouped'); hold on;

    for im = 1:numel(metricsA)
        lo = max(M(:,im) - S(:,im), 0);
        hi = min(M(:,im) + S(:,im), 1);
        errorbar(b(im).XEndPoints, M(:,im), M(:,im)-lo, hi-M(:,im), ...
            'k', 'linestyle', 'none', 'LineWidth', 1);
    end

    set(gca, 'XTick', 1:numel(groups), 'XTickLabel', groups);
    ylabel('Fraction of epochs (0–1)');
    ylim([0 1]);
    legend(labelsA, 'Location', 'best');
    title('Arousal-marked epochs');
    grid on;

    exportgraphics(gcf, fullfile(outDir, 'FIG_ArousalPct_GroupBar.png'), 'Resolution', 300);
end

% ---------- Figure 2: Visual bad percentages by stage ----------
if any(~isnan(Tsubj.VisualBadPct_WAKE)) || any(~isnan(Tsubj.VisualBadPct_N3))

    metricsV = {'VisualBadPct_WAKE','VisualBadPct_REM','VisualBadPct_N1','VisualBadPct_N2','VisualBadPct_N3'};
    labelsV  = {'WAKE','REM','N1','N2','N3'};

    M2 = nan(numel(groups), numel(metricsV));
    S2 = nan(numel(groups), numel(metricsV));

    for ig = 1:numel(groups)
        idx = (grp == groups{ig});
        for im = 1:numel(metricsV)
            x = Tsubj.(metricsV{im})(idx);
            M2(ig,im) = mean(x, 'omitnan');
            S2(ig,im) = std(x, 0, 'omitnan');
        end
    end

    figure('Color','w');
    b2 = bar(M2, 'grouped'); hold on;

    for im = 1:numel(metricsV)
        lo = max(M2(:,im) - S2(:,im), 0);
        hi = min(M2(:,im) + S2(:,im), 1);
        errorbar(b2(im).XEndPoints, M2(:,im), M2(:,im)-lo, hi-M2(:,im), ...
            'k', 'linestyle', 'none', 'LineWidth', 1);
    end

    set(gca, 'XTick', 1:numel(groups), 'XTickLabel', groups);
    ylabel('Fraction of epochs (0–1)');
    ylim([0 1]);
    legend(labelsV, 'Location', 'bestoutside');
    title('Visual-marked epochs by sleep stage');
    grid on;

    exportgraphics(gcf, fullfile(outDir, 'FIG_VisualBadPct_GroupBar.png'), 'Resolution', 300);
end

%% ===================== BASIC RECAP TABLE =====================
metrics = { ...
    'ArousalBadPct_total',  'Arousal bad (overall) %'; ...
    'ArousalBadPct_NREM',   'Arousal bad (NREM) %'; ...
    'ArousalBadPct_N3',     'Arousal bad (N3) %'; ...
    'VisualBadPct_total',   'Visual bad (overall) %'; ...
    'VisualBadPct_NREM',    'Visual bad (NREM) %'; ...
    'VisualBadPct_WAKE',    'Visual bad (WAKE) %'; ...
    'VisualBadPct_REM',     'Visual bad (REM) %'; ...
    'VisualBadPct_N1',      'Visual bad (N1) %'; ...
    'VisualBadPct_N2',      'Visual bad (N2) %'; ...
    'VisualBadPct_N3',      'Visual bad (N3) %'; ...
    'ICA_nCh_used',         'ICA channels used (n)'; ...
    'ICA_nCh_interpolated', 'Interpolated scalp channels (n)'; ...
    'ICA_nIC_bad_auto',     'ICLabel bad ICs (n)'};

nG = numel(groups);
nM = size(metrics,1);

recap = cell(nG, nM+1);
rowNames = groups(:);
colNames = cell(1, nM+1);
colNames{1} = 'N';

for im = 1:nM
    colNames{im+1} = metrics{im,2};
end

for ig = 1:nG
    idx = (grp == groups{ig});
    recap{ig,1} = sprintf('%d', sum(idx));

    for im = 1:nM
        vn = metrics{im,1};
        x  = Tsubj.(vn)(idx);
        mu = mean(x, 'omitnan');
        sd = std(x, 0, 'omitnan');

        if contains(vn, 'Pct')
            recap{ig,im+1} = sprintf('%.1f ± %.1f', 100*mu, 100*sd);
        else
            recap{ig,im+1} = sprintf('%.1f ± %.1f', mu, sd);
        end
    end
end

Trecap = cell2table(recap, 'VariableNames', matlab.lang.makeValidName(colNames));
Trecap = addvars(Trecap, string(rowNames), 'Before', 1, 'NewVariableNames', "Group");
disp(Trecap);

writetable(Trecap, fullfile(outDir, 'TABLE_QC_BasicRecap.csv'));

fprintf('Saved recap table to: %s\n', outDir);