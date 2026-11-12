%% ========================================================================
% WP2 — N3 Availability and Early–Late NREM Composition Report
%
% DESCRIPTION:
%   - Segments each subject's NREM sleep (N1+N2+N3) into 10 bins
%   - Computes per-segment counts and within-NREM composition (%N1/%N2/%N3)
%   - Summarizes Early vs Late NREM windows at the subject level
%   - Exports reviewer-friendly recap tables, including T2–T9 summaries
%   - Produces paper-style figures:
%       Fig1: group mean ± SD across segments (%N1/%N2/%N3)
%       Fig2: subject-level spaghetti plots for %N3 within each NREM segment
%       Fig3: subject-level Early vs Late barplots for %N3 within NREM
%       Fig4: N3 availability in Early vs Late sleep (counts and minutes),
%             including matched-N3 availability
%   - Saves N3 epoch indices for later matched-N3 sensitivity analyses
%
% NOTE:
%   - Internal group codes remain:
%       CTL, ADV, DYS
%   - Display labels in figures/tables are:
%       CTL, PDdys-, PDdys+
% ========================================================================

clear; clc; close all;

%% ====================== USER PARAMS =====================================
GROUP_ORDER = {'CTL','ADV','DYS'};      % internal group codes
GROUPS      = GROUP_ORDER;

NSEG        = 10;

% Early / Late definition
EARLY_SEG   = 2:5;
LATE_SEG    = 6:9;

SEG_KEEP_FOR_WP2 = 2:9;

ROOT        = '...\_processed';
SCORING_DIR = fullfile(ROOT,'infoscoringpersubject');

OUTROOT     = '...\_entropy';
OUT_DIR     = fullfile(OUTROOT,'GroupAnalysis','WP2_N3_AVAILABILITY');
ensure_dir(OUT_DIR);

% Optional: restrict subjects to those having metric files
USE_METRIC_LISTING = true;
metricSubdir  = 'SampEn_output';

PRINT_PER_SUBJECT = false;
SAVE_INDEXMAP = true;

epoch_len_s = 30;

%% ====================== COLLECT & COMPUTE ===============================
AllLong   = table();   % long format: Group/Subj/Seg + counts/percentages
AllSubj   = table();   % one row per subject: Early/Late summaries + matched counts
Missing   = {};        % log of missing/excluded subjects
IdxMap    = struct('Group',{},'Subj',{},'NSEG',{},'EARLY_SEG',{},'LATE_SEG',{}, ...
                   'idxNREM',{},'segOfNREM',{},'idxN3_Early',{},'idxN3_Late',{});

fprintf('\n=== WP2 N3 AVAILABILITY | NSEG=%d | EARLY=%s | LATE=%s ===\n', ...
    NSEG, join(string(EARLY_SEG),','), join(string(LATE_SEG),','));
fprintf('SCORING_DIR: %s\nOUT_DIR:     %s\n\n', SCORING_DIR, OUT_DIR);

for g = 1:numel(GROUPS)
    G = GROUPS{g};

    subjList = get_subjects_for_group(G, ROOT, SCORING_DIR, USE_METRIC_LISTING, metricSubdir);
    if isempty(subjList)
        fprintf('[%s] No subjects found.\n', G);
        continue;
    end

    fprintf('[%s] Subjects: %d\n', G, numel(subjList));

    for k = 1:numel(subjList)
        subj = subjList{k};

        scoringPath = fullfile(SCORING_DIR, sprintf('%s_%s_scoring.mat', G, subj));
        if ~exist(scoringPath,'file')
            Missing(end+1,:) = {G, subj, 'missing_scoring_mat'}; %#ok<AGROW>
            continue;
        end

        Tsc = load(scoringPath);
        if ~isfield(Tsc,'epochStages') || isempty(Tsc.epochStages)
            Missing(end+1,:) = {G, subj, 'missing_epochStages'}; %#ok<AGROW>
            continue;
        end

        epochStages = Tsc.epochStages(:)'; % row vector

        % ---- Per-segment table + index mapping (NREM-based) ----
        [Tseg, segInfo] = stage_counts_per_segment_NREM(epochStages, NSEG, epoch_len_s);
        if isempty(Tseg) || height(Tseg)==0
            Missing(end+1,:) = {G, subj, sprintf('NREM<NSEG (%d)',NSEG)}; %#ok<AGROW>
            continue;
        end

        % Add identifiers
        Tseg.Group = repmat(string(G), height(Tseg), 1);
        Tseg.GroupLabel = repmat(group_label_csv(G), height(Tseg), 1);
        Tseg.Subj  = repmat(string(subj), height(Tseg), 1);
        Tseg = movevars(Tseg, {'Group','GroupLabel','Subj'}, 'Before', 1);

        % ---- Early/Late summary (NREM-based) ----
        Tsum = summarize_early_late(Tseg, EARLY_SEG, LATE_SEG, epoch_len_s);
        Tsum.Group = string(G);
        Tsum.GroupLabel = group_label_csv(G);
        Tsum.Subj  = string(subj);
        Tsum = movevars(Tsum, {'Group','GroupLabel','Subj'}, 'Before', 1);

        % ---- Matched-N3 availability ----
        idxN3_E = segInfo.idxN3_in_segments(EARLY_SEG);
        idxN3_L = segInfo.idxN3_in_segments(LATE_SEG);
        nN3_E   = numel(idxN3_E);
        nN3_L   = numel(idxN3_L);
        nMatch  = min(nN3_E, nN3_L);

        Tsum.nN3_E = nN3_E;
        Tsum.nN3_L = nN3_L;
        Tsum.nN3_match = nMatch;
        Tsum.frac_match_vs_E = safe_div(nMatch, nN3_E);
        Tsum.frac_match_vs_L = safe_div(nMatch, nN3_L);
        Tsum.minN3_E     = nN3_E * epoch_len_s / 60;
        Tsum.minN3_L     = nN3_L * epoch_len_s / 60;
        Tsum.minN3_match = nMatch * epoch_len_s / 60;

        % Accumulate
        AllLong = [AllLong; Tseg]; %#ok<AGROW>
        AllSubj = [AllSubj; Tsum]; %#ok<AGROW>

        % Index map for later matched-N3 on SampEn
        if SAVE_INDEXMAP
            s0.Group = G;
            s0.Subj  = subj;
            s0.NSEG  = NSEG;
            s0.EARLY_SEG = EARLY_SEG;
            s0.LATE_SEG  = LATE_SEG;
            s0.idxNREM   = segInfo.idxNREM;
            s0.segOfNREM = segInfo.segOfNREM;
            s0.idxN3_Early = idxN3_E;
            s0.idxN3_Late  = idxN3_L;
            IdxMap(end+1) = s0; %#ok<AGROW>
        end

        if PRINT_PER_SUBJECT
            print_segment_report(subj, Tseg, EARLY_SEG, LATE_SEG);
        end
    end
end

%% ====================== POST-PROCESS / ORDERING =========================
if ~isempty(AllLong)
    AllLong.Group = categorical(AllLong.Group, GROUP_ORDER, 'Ordinal', true);
end
if ~isempty(AllSubj)
    AllSubj.Group = categorical(AllSubj.Group, GROUP_ORDER, 'Ordinal', true);
end

%% ====================== SAVE OUTPUTS ====================================
if ~isempty(AllLong)
    out1 = fullfile(OUT_DIR, sprintf('WP2_SegmentCounts_Long_NSEG=%d.csv', NSEG));
    writetable(AllLong, out1);
    fprintf('\nSaved: %s\n', out1);
else
    fprintf('\nNo valid subjects -> AllLong empty.\n');
end

if ~isempty(AllSubj)
    out2 = fullfile(OUT_DIR, sprintf('WP2_SegmentCounts_EarlyLate_NSEG=%d.csv', NSEG));
    writetable(AllSubj, out2);
    fprintf('Saved: %s\n', out2);
end

% Group summary: mean ± SD across subjects per group × segment
if ~isempty(AllLong)
    Gsum = make_group_summary(AllLong);
    out3 = fullfile(OUT_DIR, sprintf('WP2_SegmentCounts_GroupSummary_NSEG=%d.csv', NSEG));
    writetable(Gsum, out3);
    fprintf('Saved: %s\n', out3);

    % Reviewer-friendly: only T2–T9
    Gsum_T2T9 = Gsum(ismember(Gsum.Seg, SEG_KEEP_FOR_WP2), :);
    out3b = fullfile(OUT_DIR, sprintf('WP2_GroupSummary_T2T9_NSEG=%d.csv', NSEG));
    writetable(Gsum_T2T9, out3b);
    fprintf('Saved: %s\n', out3b);
end

% Missing log
if ~isempty(Missing)
    Tmiss = cell2table(Missing, 'VariableNames', {'Group','Subj','Reason'});
    Tmiss.GroupLabel = strings(height(Tmiss),1);
    for ii = 1:height(Tmiss)
        Tmiss.GroupLabel(ii) = group_label_csv(Tmiss.Group{ii});
    end
    Tmiss = movevars(Tmiss, 'GroupLabel', 'After', 'Group');

    outM  = fullfile(OUT_DIR, sprintf('WP2_MissingOrExcluded_NSEG=%d.csv', NSEG));
    writetable(Tmiss, outM);
    fprintf('Saved: %s\n', outM);
end

% Index map MAT
if SAVE_INDEXMAP && ~isempty(IdxMap)
    outMat = fullfile(OUT_DIR, sprintf('WP2_N3_IndexMap_NSEG=%d.mat', NSEG));
    save(outMat, 'IdxMap');
    fprintf('Saved: %s\n', outMat);
end

fprintf('\n✅ WP2 tables saved.\n');

%% ====================== PLOTS (PAPER-LIKE) ==============================
if ~isempty(AllLong)

    % Fig1: group mean ± SD across segments (N1/N2/N3 %)
    plot_groups_mean_sd(AllLong, OUT_DIR, NSEG);

    % Fig2: subject-level spaghetti plots — %N3 within NREM segment
    plot_subjects_spaghetti_colored(AllLong, OUT_DIR, NSEG, 'N3', GROUP_ORDER);

    % Fig3: subject-level Early vs Late bars — %N3 within NREM
    plot_subjects_early_late_bars(AllLong, OUT_DIR, NSEG, 'N3', EARLY_SEG, LATE_SEG, GROUP_ORDER);
end

if ~isempty(AllSubj)
    % Fig4: N3 availability Early vs Late (counts + minutes + matched)
    plot_n3_availability_counts(AllSubj, OUT_DIR, GROUP_ORDER);
end

fprintf('\n✅ Done.\n');

%% ===================== BUILD TABLES FROM AllSubj ========================
% Expected columns in AllSubj:
% Group, Subj, pN3_E, pN3_L, etc.

GROUP_ORDER = {'CTL','ADV','DYS'};
AllSubj.Group = categorical(AllSubj.Group, GROUP_ORDER, 'Ordinal', true);

%% ===================== 1) GROUP TABLE (mean ± SD) =======================
Tgroup = table();
for gi = 1:numel(GROUP_ORDER)
    G = GROUP_ORDER{gi};
    Tg = AllSubj(AllSubj.Group==G, :);
    if isempty(Tg), continue; end

    row = table();
    row.Group = string(group_label_csv(G));
    row.nSubj = height(Tg);

    row.mean_N3early = mean(Tg.pN3_E, 'omitnan');
    row.sd_N3early   = std(Tg.pN3_E, 0, 'omitnan');

    row.mean_N3late  = mean(Tg.pN3_L, 'omitnan');
    row.sd_N3late    = std(Tg.pN3_L, 0, 'omitnan');

    d = Tg.pN3_L - Tg.pN3_E;
    row.mean_Delta   = mean(d, 'omitnan');
    row.sd_Delta     = std(d, 0, 'omitnan');

    Tgroup = [Tgroup; row]; %#ok<AGROW>
end

Tgroup_paper = Tgroup(:,{'Group','nSubj','mean_N3early','sd_N3early','mean_N3late','sd_N3late'});
Tgroup_paper.N3early_meanSD = compose('%.1f ± %.1f', Tgroup.mean_N3early, Tgroup.sd_N3early);
Tgroup_paper.N3late_meanSD  = compose('%.1f ± %.1f', Tgroup.mean_N3late,  Tgroup.sd_N3late);
Tgroup_paper = movevars(Tgroup_paper, {'N3early_meanSD','N3late_meanSD'}, 'After', 'nSubj');
Tgroup_paper(:,{'mean_N3early','sd_N3early','mean_N3late','sd_N3late'}) = [];

disp('=== GROUP SUMMARY (%N3 within NREM) ===');
disp(Tgroup_paper);

out_group = fullfile(OUT_DIR, 'TABLE_GroupSummary_N3earlylate.csv');
writetable(Tgroup_paper, out_group);

%% ===================== 2) PER-SUBJECT TABLES ============================
make_subj_table = @(G) local_make_subject_table(AllSubj, G);

T_CTL = make_subj_table('CTL');
T_ADV = make_subj_table('ADV');
T_DYS = make_subj_table('DYS');

disp('=== CTL per-subject ===');      disp(T_CTL);
disp('=== PDdys- per-subject ===');   disp(T_ADV);
disp('=== PDdys+ per-subject ===');   disp(T_DYS);

writetable(T_CTL, fullfile(OUT_DIR,'TABLE_Subjects_CTL_N3earlylate.csv'));
writetable(T_ADV, fullfile(OUT_DIR,'TABLE_Subjects_PDdysMinus_N3earlylate.csv'));
writetable(T_DYS, fullfile(OUT_DIR,'TABLE_Subjects_PDdysPlus_N3earlylate.csv'));

fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', ...
    out_group, ...
    fullfile(OUT_DIR,'TABLE_Subjects_CTL_N3earlylate.csv'), ...
    fullfile(OUT_DIR,'TABLE_Subjects_PDdysMinus_N3earlylate.csv'), ...
    fullfile(OUT_DIR,'TABLE_Subjects_PDdysPlus_N3earlylate.csv'));

%% ===================== local function ===================================
function Tout = local_make_subject_table(AllSubj, G)
Tg = AllSubj(AllSubj.Group==G, :);
if isempty(Tg)
    Tout = table();
    return;
end

Tout = table();
Tout.Group     = repmat(group_label_csv(G), height(Tg), 1);
Tout.Subj      = string(Tg.Subj);
Tout.N3early   = Tg.pN3_E;
Tout.N3late    = Tg.pN3_L;
Tout.Delta_LmE = Tout.N3late - Tout.N3early;

[~,ord] = sort(Tout.Delta_LmE, 'descend', 'MissingPlacement','last');
Tout = Tout(ord,:);
end

%% ========================================================================
%%                                FUNCTIONS
%% ========================================================================

function y = safe_div(a,b)
if b==0
    y = NaN;
else
    y = a./b;
end
end

function ensure_dir(d)
if ~exist(d,'dir')
    mkdir(d);
end
end

function lbl = group_label_csv(G)
% Plain-text label for tables / CSV / axis tick labels
G = upper(string(G));
switch G
    case "ADV"
        lbl = "PDdys-";
    case "DYS"
        lbl = "PDdys+";
    otherwise
        lbl = "CTL";
end
end

function lbl = group_label_tex(G)
% TeX label for plot titles / legends
G = upper(string(G));
switch G
    case "ADV"
        lbl = 'PDdys^{-}';
    case "DYS"
        lbl = 'PDdys^{+}';
    otherwise
        lbl = 'CTL';
end
end

function subjList = get_subjects_for_group(G, ROOT, SCORING_DIR, USE_METRIC_LISTING, metricSubdir)
% Prefer metric listing if requested; otherwise use scoring files.
subjList = {};

if USE_METRIC_LISTING
    p = fullfile(ROOT, metricSubdir, G);
    if exist(p,'dir')
        d = dir(fullfile(p, '*_eeg_minimal_*.mat'));
        for i = 1:numel(d)
            subjList{end+1} = extract_subj_id(d(i).name); %#ok<AGROW>
        end
        subjList = unique(subjList, 'stable');
        return;
    end
end

d = dir(fullfile(SCORING_DIR, sprintf('%s_*_scoring.mat', G)));
for i = 1:numel(d)
    subjList{end+1} = extract_subj_id_from_scoring(d(i).name, G); %#ok<AGROW>
end
subjList = unique(subjList, 'stable');
end

function subj = extract_subj_id(fname_metric)
tok = regexp(fname_metric,'^(.*?)_eeg_minimal','tokens','once');
if isempty(tok)
    subj = erase(fname_metric,'.mat');
else
    subj = tok{1};
end
end

function subj = extract_subj_id_from_scoring(fname_scoring, G)
pat = sprintf('^%s_(.*)_scoring\\.mat$', regexptranslate('escape',G));
tok = regexp(fname_scoring, pat, 'tokens','once');
if isempty(tok)
    subj = erase(fname_scoring, '_scoring.mat');
else
    subj = tok{1};
end
end

function [Tseg, segInfo] = stage_counts_per_segment_NREM(epochStages, NSEG, epoch_len_s)
% Segment NREM epochs (N1+N2+N3) into NSEG bins.
% Count n and percentages within each bin.
% Add minutes using epoch_len_s.

epochStages = epochStages(:)';  % row
idxNREM = find(epochStages==1 | epochStages==2 | epochStages==3);

if numel(idxNREM) < NSEG
    Tseg = table();
    segInfo = struct();
    return;
end

edges = round(linspace(1, numel(idxNREM)+1, NSEG+1));

Seg  = (1:NSEG)';
nTot = zeros(NSEG,1);
nN1  = zeros(NSEG,1);
nN2  = zeros(NSEG,1);
nN3  = zeros(NSEG,1);

segOfNREM = zeros(1, numel(idxNREM));
idxN3_bySeg = cell(1,NSEG);

for s = 1:NSEG
    range = edges(s):edges(s+1)-1;
    segNREMpos = range;
    segIdx = idxNREM(segNREMpos);
    st     = epochStages(segIdx);

    nTot(s) = numel(segIdx);
    nN1(s)  = sum(st==1);
    nN2(s)  = sum(st==2);
    nN3(s)  = sum(st==3);

    segOfNREM(segNREMpos) = s;
    idxN3_bySeg{s} = segIdx(st==3);
end

pN1 = 100 * nN1 ./ max(nTot,1);
pN2 = 100 * nN2 ./ max(nTot,1);
pN3 = 100 * nN3 ./ max(nTot,1);

minTot = nTot * epoch_len_s / 60;
minN1  = nN1  * epoch_len_s / 60;
minN2  = nN2  * epoch_len_s / 60;
minN3  = nN3  * epoch_len_s / 60;

Tseg = table(Seg, nTot, nN1, nN2, nN3, pN1, pN2, pN3, ...
             minTot, minN1, minN2, minN3, ...
    'VariableNames', {'Seg','nTot','nN1','nN2','nN3','pN1','pN2','pN3', ...
                      'minTot','minN1','minN2','minN3'});

segInfo.idxNREM = idxNREM;
segInfo.segOfNREM = segOfNREM;
segInfo.idxN3_bySeg = idxN3_bySeg;
segInfo.idxN3_in_segments = @(segList) cat(2, idxN3_bySeg{segList});
end

function Tsum = summarize_early_late(Tseg, earlySeg, lateSeg, epoch_len_s)
% Subject-level summary (EARLY vs LATE totals).
% Percentages are within NREM.
% Adds minutes.

E = ismember(Tseg.Seg, earlySeg);
L = ismember(Tseg.Seg, lateSeg);

nE   = sum(Tseg.nTot(E));
nL   = sum(Tseg.nTot(L));

nE1  = sum(Tseg.nN1(E)); nL1 = sum(Tseg.nN1(L));
nE2  = sum(Tseg.nN2(E)); nL2 = sum(Tseg.nN2(L));
nE3  = sum(Tseg.nN3(E)); nL3 = sum(Tseg.nN3(L));

pE1  = 100*nE1/max(nE,1); pL1 = 100*nL1/max(nL,1);
pE2  = 100*nE2/max(nE,1); pL2 = 100*nL2/max(nL,1);
pE3  = 100*nE3/max(nE,1); pL3 = 100*nL3/max(nL,1);

minTot_E = nE  * epoch_len_s / 60;
minN3_E  = nE3 * epoch_len_s / 60;
minTot_L = nL  * epoch_len_s / 60;
minN3_L  = nL3 * epoch_len_s / 60;

Tsum = table(nE,nE1,nE2,nE3,pE1,pE2,pE3, ...
             nL,nL1,nL2,nL3,pL1,pL2,pL3, ...
             minTot_E,minN3_E,minTot_L,minN3_L, ...
    'VariableNames', {'nTot_E','nN1_E','nN2_E','nN3_E','pN1_E','pN2_E','pN3_E', ...
                      'nTot_L','nN1_L','nN2_L','nN3_L','pN1_L','pN2_L','pN3_L', ...
                      'minTot_E','minN3_E','minTot_L','minN3_L'});
end

function print_segment_report(subj, Tseg, earlySeg, lateSeg)
fprintf('  %s | per-segment (NREM-based):\n', subj);
for i = 1:height(Tseg)
    fprintf('    Seg %02d: n=%3d | N1=%3d (%.1f%%)  N2=%3d (%.1f%%)  N3=%3d (%.1f%%)\n', ...
        Tseg.Seg(i), Tseg.nTot(i), ...
        Tseg.nN1(i), Tseg.pN1(i), ...
        Tseg.nN2(i), Tseg.pN2(i), ...
        Tseg.nN3(i), Tseg.pN3(i));
end
E = ismember(Tseg.Seg, earlySeg);
L = ismember(Tseg.Seg, lateSeg);
nE = sum(Tseg.nTot(E)); nL = sum(Tseg.nTot(L));
nE3 = sum(Tseg.nN3(E)); nL3 = sum(Tseg.nN3(L));
fprintf('    EARLY seg (%s): N3=%d/%d (%.1f%% of NREM)\n', join(string(earlySeg),','), nE3, nE, 100*nE3/max(nE,1));
fprintf('    LATE  seg (%s): N3=%d/%d (%.1f%% of NREM)\n', join(string(lateSeg),','),  nL3, nL, 100*nL3/max(nL,1));
end

function Gsum = make_group_summary(AllLong)
% Mean ± SD across subjects per Group × Seg for pN1/pN2/pN3 and counts.

T = AllLong(:, {'Group','GroupLabel','Subj','Seg','nTot','nN1','nN2','nN3','pN1','pN2','pN3'});
groups = categories(T.Group);
segs   = unique(T.Seg);
Gsum   = table();

for gi = 1:numel(groups)
    G = groups{gi};
    for si = 1:numel(segs)
        s = segs(si);

        Ts = T(T.Group==G & T.Seg==s, :);
        if isempty(Ts), continue; end

        nSubj = numel(unique(Ts.Subj));
        row = table();
        row.Group = string(G);
        row.GroupLabel = group_label_csv(G);
        row.Seg   = s;
        row.nSubj = nSubj;

        row.mean_pN1 = mean(Ts.pN1,'omitnan'); row.sd_pN1 = std(Ts.pN1,0,'omitnan');
        row.mean_pN2 = mean(Ts.pN2,'omitnan'); row.sd_pN2 = std(Ts.pN2,0,'omitnan');
        row.mean_pN3 = mean(Ts.pN3,'omitnan'); row.sd_pN3 = std(Ts.pN3,0,'omitnan');

        row.mean_nTot = mean(Ts.nTot,'omitnan'); row.sd_nTot = std(Ts.nTot,0,'omitnan');
        row.mean_nN3  = mean(Ts.nN3,'omitnan');  row.sd_nN3  = std(Ts.nN3,0,'omitnan');

        Gsum = [Gsum; row]; %#ok<AGROW>
    end
end
end

%% =============================== PLOTS ==================================

function plot_groups_mean_sd(AllLong, OUT_DIR, NSEG)
% One figure, 3 panels: pN1, pN2, pN3. Lines = groups, band = ± SD.

stages = {'pN1','pN2','pN3'};
stageTitles = {'% N1 within NREM','% N2 within NREM','% N3 within NREM'};

groups = categories(AllLong.Group);
cols = lines(numel(groups));

fig = figure('Color','w','Position',[80 80 1600 480]);
t = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for si = 1:3
    ax = nexttile(t,si); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    for gi = 1:numel(groups)
        G = groups{gi};
        Tg = AllLong(AllLong.Group==G, :);

        mu = nan(1,NSEG);
        sd = nan(1,NSEG);
        for s = 1:NSEG
            y = Tg.(stages{si})(Tg.Seg==s);
            mu(s) = mean(y,'omitnan');
            sd(s) = std(y,0,'omitnan');
        end

        x = 1:NSEG;
        good = ~isnan(mu) & ~isnan(sd);
        if any(good)
            xx = x(good);
            patch(ax, [xx fliplr(xx)], [(mu(good)-sd(good)) fliplr(mu(good)+sd(good))], ...
                cols(gi,:), 'FaceAlpha',0.15,'EdgeColor','none');
        end

        plot(ax, x, mu, '-o', ...
            'LineWidth', 2, ...
            'MarkerSize', 4, ...
            'DisplayName', group_label_tex(G), ...
            'Color', cols(gi,:));
    end

    xlim(ax,[1 NSEG]); ylim(ax,[0 100]);
    xlabel(ax,'Segment (NREM-based)');
    ylabel(ax,stageTitles{si});
    title(ax,stageTitles{si},'FontWeight','bold');
end

legend(nexttile(t,1), 'Location','northoutside','Orientation','horizontal','Box','off');
save_fig(fig, OUT_DIR, sprintf('FIG1_GroupMeanSD_pN1pN2pN3_NSEG=%d', NSEG));
end

function plot_subjects_spaghetti_colored(AllLong, OUT_DIR, NSEG, WHICH_STAGE, GROUP_ORDER)
% Rows = groups (fixed order). Each subject has its own color. Group mean = thick black.

if nargin<5
    GROUP_ORDER = categories(AllLong.Group);
end

switch upper(WHICH_STAGE)
    case 'N1'
        varY = 'pN1'; ttl = '% N1 within NREM';
    case 'N2'
        varY = 'pN2'; ttl = '% N2 within NREM';
    otherwise
        varY = 'pN3'; ttl = '% N3 within NREM';
end

groups = GROUP_ORDER;
figH = 280*numel(groups);
fig  = figure('Color','w','Position',[80 80 1400 figH]);
t    = tiledlayout(numel(groups),1,'Padding','compact','TileSpacing','compact');

for gi = 1:numel(groups)
    G = groups{gi};
    ax = nexttile(t,gi); hold(ax,'on'); grid(ax,'on'); box(ax,'on');

    Tg   = AllLong(AllLong.Group==G,:);
    subs = unique(Tg.Subj);
    nS   = numel(subs);
    cmap = lines(max(nS,1));

    for k = 1:nS
        S  = subs(k);
        Ts = Tg(Tg.Subj==S, :);

        y = nan(1,NSEG);
        for s = 1:NSEG
            row = Ts(Ts.Seg==s,:);
            if ~isempty(row)
                y(s) = row.(varY)(1);
            end
        end
        plot(ax, 1:NSEG, y, '-', 'LineWidth',1.4, 'Color', cmap(k,:));
    end

    mu = nan(1,NSEG);
    for s = 1:NSEG
        ys = Tg.(varY)(Tg.Seg==s);
        mu(s) = mean(ys,'omitnan');
    end
    plot(ax, 1:NSEG, mu, '-o', 'LineWidth',3.5, 'MarkerSize',4, 'Color', [0 0 0]);

    xlim(ax,[1 NSEG]); ylim(ax,[0 100]);
    ylabel(ax, group_label_tex(G), 'FontWeight','bold', 'Interpreter','tex');
    if gi==numel(groups)
        xlabel(ax,'Segment (NREM-based)');
    end
end

title(t, sprintf('FIG2 Subjects spaghetti (colored) — %s per segment', ttl), 'FontWeight','bold');
save_fig(fig, OUT_DIR, sprintf('FIG2_SubjectsSpaghettiColored_%s_NSEG=%d', upper(WHICH_STAGE), NSEG));
end

function plot_n3_availability_counts(AllSubj, OUT_DIR, GROUP_ORDER)
% FIG4: N3 availability Early vs Late (counts and minutes) + matched info (2x2)

if nargin<3
    GROUP_ORDER = categories(AllSubj.Group);
end

groups = GROUP_ORDER;
cols = lines(numel(groups));
xlab = cellfun(@(x) char(group_label_csv(x)), groups, 'UniformOutput', false);

fig = figure('Color','w','Position',[80 80 1400 800]);
t = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% ---------- A) Counts paired ----------
ax1 = nexttile(t,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
for gi=1:numel(groups)
    G = groups{gi};
    Tg = AllSubj(AllSubj.Group==G,:);
    x1 = gi-0.12; x2 = gi+0.12;
    for i=1:height(Tg)
        plot(ax1,[x1 x2],[Tg.nN3_E(i) Tg.nN3_L(i)],'-','Color',cols(gi,:), 'LineWidth',1.2);
        plot(ax1,x1,Tg.nN3_E(i),'o','Color',cols(gi,:), 'MarkerFaceColor',cols(gi,:));
        plot(ax1,x2,Tg.nN3_L(i),'o','Color',cols(gi,:), 'MarkerFaceColor',cols(gi,:));
    end
end
set(ax1,'XTick',1:numel(groups),'XTickLabel',xlab);
ylabel(ax1,'N3 epochs (count)');
title(ax1,'FIG4A N3 availability — Early vs Late (epochs)','FontWeight','bold');

% ---------- B) Minutes paired ----------
ax2 = nexttile(t,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
for gi=1:numel(groups)
    G = groups{gi};
    Tg = AllSubj(AllSubj.Group==G,:);
    x1 = gi-0.12; x2 = gi+0.12;
    for i=1:height(Tg)
        plot(ax2,[x1 x2],[Tg.minN3_E(i) Tg.minN3_L(i)],'-','Color',cols(gi,:), 'LineWidth',1.2);
        plot(ax2,x1,Tg.minN3_E(i),'o','Color',cols(gi,:), 'MarkerFaceColor',cols(gi,:));
        plot(ax2,x2,Tg.minN3_L(i),'o','Color',cols(gi,:), 'MarkerFaceColor',cols(gi,:));
    end
end
set(ax2,'XTick',1:numel(groups),'XTickLabel',xlab);
ylabel(ax2,'N3 minutes');
title(ax2,'FIG4B N3 availability — Early vs Late (minutes)','FontWeight','bold');

% ---------- C) Matched counts ----------
ax3 = nexttile(t,3); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
mu = nan(1,numel(groups)); sd = nan(1,numel(groups));
for gi=1:numel(groups)
    G = groups{gi};
    Tg = AllSubj(AllSubj.Group==G,:);
    mu(gi) = mean(Tg.nN3_match,'omitnan');
    sd(gi) = std(Tg.nN3_match,0,'omitnan');
    jitter = (rand(height(Tg),1)-0.5)*0.18;
    scatter(ax3, gi + jitter, Tg.nN3_match, 36, 'filled');
end
errorbar(ax3, 1:numel(groups), mu, sd, 'k', 'LineStyle','none','LineWidth',1.5);
set(ax3,'XTick',1:numel(groups),'XTickLabel',xlab);
ylabel(ax3,'Matched N3 epochs = min(Early, Late)');
title(ax3,'FIG4C Matched-N3 availability (epochs)','FontWeight','bold');

% ---------- D) Matched minutes ----------
ax4 = nexttile(t,4); hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');
mu = nan(1,numel(groups)); sd = nan(1,numel(groups));
for gi=1:numel(groups)
    G = groups{gi};
    Tg = AllSubj(AllSubj.Group==G,:);
    mu(gi) = mean(Tg.minN3_match,'omitnan');
    sd(gi) = std(Tg.minN3_match,0,'omitnan');
    jitter = (rand(height(Tg),1)-0.5)*0.18;
    scatter(ax4, gi + jitter, Tg.minN3_match, 36, 'filled');
end
errorbar(ax4, 1:numel(groups), mu, sd, 'k', 'LineStyle','none','LineWidth',1.5);
set(ax4,'XTick',1:numel(groups),'XTickLabel',xlab);
ylabel(ax4,'Matched N3 minutes = min(Early, Late)');
title(ax4,'FIG4D Matched-N3 availability (minutes)','FontWeight','bold');

save_fig(fig, OUT_DIR, 'FIG4_N3Availability_EarlyLate_epochs_minutes_and_Matched');
end

function save_fig(fig, outdir, basename)
% Save PNG (300 dpi) + vector PDF
png = fullfile(outdir, [basename '.png']);
pdf = fullfile(outdir, [basename '.pdf']);
exportgraphics(fig, png, 'Resolution', 300);
exportgraphics(fig, pdf, 'ContentType','vector');
fprintf('Saved: %s\n', png);
fprintf('Saved: %s\n', pdf);
end

function plot_subjects_early_late_bars(AllLong, OUT_DIR, NSEG, WHICH_STAGE, earlySeg, lateSeg, GROUP_ORDER)
% FIG3
% - one subplot per group (vertical layout)
% - paired bars Early vs Late per subject
% - subjects sorted by delta (Late - Early)
% - compact spacing between subjects
% - group labels: CTL, PDdys^{-}, PDdys^{+}
% - subject labels on x-axis: PDxxx -> IDxxx
% - mean ± SD box at top-right of each subplot
% - paper-like style
% - export at 600 dpi

if nargin < 7 || isempty(GROUP_ORDER)
    GROUP_ORDER = categories(AllLong.Group);
end

switch upper(WHICH_STAGE)
    case 'N1'
        varY = 'pN1';
        ttl  = '%N1 within NREM';
    case 'N2'
        varY = 'pN2';
        ttl  = '%N2 within NREM';
    otherwise
        varY = 'pN3';
        ttl  = '%N3 within NREM';
end

groups = GROUP_ORDER;

% ---------- Global y-limit ----------
allVals = [];
for gi = 1:numel(groups)
    G  = groups{gi};
    Tg = AllLong(AllLong.Group==G,:);
    if isempty(Tg), continue; end

    subs = unique(Tg.Subj);
    E = nan(numel(subs),1);
    L = nan(numel(subs),1);

    for k = 1:numel(subs)
        S  = subs(k);
        Ts = Tg(Tg.Subj==S,:);

        y = nan(1,NSEG);
        for s = 1:NSEG
            row = Ts(Ts.Seg==s,:);
            if ~isempty(row)
                y(s) = row.(varY)(1);
            end
        end

        E(k) = mean(y(earlySeg),'omitnan');
        L(k) = mean(y(lateSeg),'omitnan');
    end

    allVals = [allVals; E(:); L(:)]; %#ok<AGROW>
end

allVals = allVals(isfinite(allVals));
if isempty(allVals)
    warning('No valid values to plot.');
    return;
end

yMax = ceil(max(allVals)/10)*10;
yMax = max(yMax, 80);
yMax = min(yMax, 100);

% ---------- Figure ----------
figH = 360 * numel(groups);
figW = 1250;
fig  = figure('Color','w','Position',[80 80 figW figH]);
t    = tiledlayout(numel(groups),1,'Padding','compact','TileSpacing','compact');

for gi = 1:numel(groups)
    G = groups{gi};
    Gdisp = group_label_tex(G);

    Tg = AllLong(AllLong.Group==G,:);
    ax = nexttile(t,gi);

    if isempty(Tg)
        axis(ax,'off');
        continue;
    end

    subs = unique(Tg.Subj);

    E = nan(numel(subs),1);
    L = nan(numel(subs),1);

    for k = 1:numel(subs)
        S  = subs(k);
        Ts = Tg(Tg.Subj==S,:);

        y = nan(1,NSEG);
        for s = 1:NSEG
            row = Ts(Ts.Seg==s,:);
            if ~isempty(row)
                y(s) = row.(varY)(1);
            end
        end

        E(k) = mean(y(earlySeg),'omitnan');
        L(k) = mean(y(lateSeg),'omitnan');
    end

    % Sort by delta (Late - Early)
    d = L - E;
    [~,ord] = sort(d,'descend','MissingPlacement','last');
    subs = subs(ord);
    E    = E(ord);
    L    = L(ord);

    % Convert subject labels: PD025 -> ID025
    subsDisp = regexprep(string(subs), '^[Pp][Dd]', 'ID');

    % Group summary
    muE = mean(E,'omitnan');
    sdE = std(E,0,'omitnan');
    muL = mean(L,'omitnan');
    sdL = std(L,0,'omitnan');

    hold(ax,'on');

    % Compact x positions
    step  = 0.62;
    X     = 1 + (0:numel(subs)-1) * step;

    delta = 0.085;
    bw    = 0.13;

    % Paper-like styling
    ax.Color      = 'w';
    ax.Box        = 'off';
    ax.LineWidth  = 1.2;
    ax.FontSize   = 16;
    ax.FontWeight = 'bold';
    ax.TickDir    = 'out';
    ax.TickLength = [0.010 0.010];
    ax.XColor     = [0.15 0.15 0.15];
    ax.YColor     = [0.15 0.15 0.15];
    ax.Layer      = 'top';
    ax.XGrid      = 'off';
    ax.YGrid      = 'on';
    ax.GridAlpha  = 0.18;
    ax.GridColor  = [0 0 0];
    ax.YMinorGrid = 'off';

    b1 = bar(ax, X-delta, E, bw, ...
        'FaceColor',[0 0.4470 0.7410], ...
        'EdgeColor',[0.20 0.20 0.20], ...
        'LineWidth',0.9);

    b2 = bar(ax, X+delta, L, bw, ...
        'FaceColor',[0.8500 0.3250 0.0980], ...
        'EdgeColor',[0.20 0.20 0.20], ...
        'LineWidth',0.9);

    xlim(ax,[X(1)-0.28, X(end)+0.28]);
    ylim(ax,[0 yMax]);

    ylabel(ax, ttl, 'FontWeight','bold', 'FontSize',18);

    title(ax, sprintf('%s — Early vs Late (per subject)', Gdisp), ...
        'FontWeight','bold', ...
        'FontSize',19, ...
        'Interpreter','tex');

    xtlbl = cellstr(subsDisp);
    set(ax, 'XTick', X, ...
            'XTickLabel', xtlbl, ...
            'FontSize', 15, ...
            'FontWeight', 'bold');
    xtickangle(ax,45);

    ax.XAxis.FontSize   = 15;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontSize   = 15;
    ax.YAxis.FontWeight = 'bold';

    if gi == 1
        lg = legend(ax, [b1 b2], {'Early','Late'}, ...
            'Location','northoutside', ...
            'Orientation','horizontal', ...
            'Box','off');
        lg.FontSize   = 16;
        lg.FontWeight = 'bold';
    end

    txt = sprintf('Early: %.1f \\pm %.1f\nLate: %.1f \\pm %.1f', ...
        muE, sdE, muL, sdL);

    text(ax, 0.985, 0.90, txt, ...
        'Units','normalized', ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'FontWeight','bold', ...
        'FontSize',13.5, ...
        'BackgroundColor','w', ...
        'EdgeColor',[0.55 0.55 0.55], ...
        'LineWidth',0.8, ...
        'Margin',8);
end

sgtitle(t, ttl, 'FontWeight','bold', 'FontSize',22);

basename = sprintf('FIG3_SubjectsEarlyLateBars_%s_NSEG=%d', upper(WHICH_STAGE), NSEG);
pngFile  = fullfile(OUT_DIR, [basename '.png']);
pdfFile  = fullfile(OUT_DIR, [basename '.pdf']);

exportgraphics(fig, pngFile, 'Resolution', 600);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

fprintf('Saved: %s\n', pngFile);
fprintf('Saved: %s\n', pdfFile);
end