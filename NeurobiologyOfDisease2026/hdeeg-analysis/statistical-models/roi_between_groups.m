%% ========================================================================
% RAW EARLY/LATE BETWEEN-GROUP ANALYSIS — SWA CLASSIC + SampEn
%
% This script performs between-group analyses on raw Early and Late values
% for two metrics:
%   - SWA classic
%   - SampEn
%
% Expected input tables must contain:
%   - Panel
%   - Group
%   - Id
%   - Early
%   - Late
%
% Statistical workflow for each metric, panel, and phase:
%   1) Group descriptives
%   2) Normality within each group
%   3) Variance homogeneity (Levene)
%   4) Residual normality from LM(Outcome ~ Group)
%   5) Omnibus test:
%        - one-way ANOVA if assumptions are satisfied
%        - Kruskal-Wallis otherwise
%   6) Pairwise post-hoc:
%        - after ANOVA: ttest2 + CI95 + Hedges g
%        - after KW   : ranksum + Hodges-Lehmann + CI95 + rank-biserial
%   7) Linear model check: fitlm(Outcome ~ Group)
%
% Holm correction is applied within each panel × phase across:
%   CTL vs ADV, CTL vs DYS, ADV vs DYS
%
% Outputs:
%   - REPORT_DESCRIPTIVES.csv
%   - REPORT_OMNIBUS.csv
%   - REPORT_PAIRWISE.csv
%   - REPORT_LM.csv
%   - METHODS_RAW_EARLY_LATE_GROUPONLY.txt
%   - combined figures (.png, .pdf)
%% ========================================================================
clear; clc; close all; rng(42);

%% ===================== USER SETTINGS ====================================
F_SWA = '...\SWA_CLASSIC_EARLY_LATE.mat';
F_SE  = '...\SampEn_EARLY_LATE.mat';

OUTDIR = '...\REPORTS_RAW_EARLY_LATE_GROUPONLY';
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

% Optional path to swtest.m. If unavailable, lillietest is used.
SWTEST_MFILE = '...\swtest\swtest.m';

GROUPS        = {'CTL','ADV','DYS'};
GROUP_DISPLAY = {'CTL','PDdys^{-}','PDdys^{+}'};

PANELS        = {'OVERALL','FRONTAL','CENTPAR','OCCIP'};
PANEL_DISPLAY = {'Overall','Frontal','Centro-parietal','Occipital'};

PHASE_VARS    = {'Early','Late'};
PHASE_DISPLAY = {'Early','Late'};

ALPHA = 0.05;

% If true, ANOVA is used only when:
% all groups normal + equal variances + residual normality + enough N
STRICT_PARAMETRIC_RULE = true;

% Holm correction within each panel × phase across 3 pairwise comparisons
NcomparisonsHolm = 3;

%% ===================== PATHS / DEBUG ====================================
fprintf('\n============================================================\n');
fprintf(' RAW EARLY/LATE BETWEEN-GROUP ANALYSIS START\n');
fprintf('============================================================\n');

if ~isempty(SWTEST_MFILE) && exist(SWTEST_MFILE,'file')
    addpath(fileparts(SWTEST_MFILE));
    fprintf('swtest added from:\n%s\n', SWTEST_MFILE);
else
    if exist('swtest','file') == 2
        fprintf('swtest already available on MATLAB path.\n');
    else
        fprintf('swtest not found -> fallback to lillietest.\n');
    end
end

fprintf('Output folder:\n%s\n', OUTDIR);
fprintf('Holm setting: NcomparisonsHolm = %d\n', NcomparisonsHolm);

%% ===================== LOAD TABLES ======================================
T_SWA = load_table_from_saved_mat(F_SWA);
T_SE  = load_table_from_saved_mat(F_SE);

T_SWA.Group = categorical(T_SWA.Group, GROUPS, GROUPS);
T_SE.Group  = categorical(T_SE.Group,  GROUPS, GROUPS);

%% ===================== RUN ANALYSIS =====================================
[Desc_SWA, Omni_SWA, Pair_SWA, LM_SWA] = run_raw_value_analysis( ...
    T_SWA, 'SWA_classic', GROUPS, GROUP_DISPLAY, ...
    PANELS, PANEL_DISPLAY, PHASE_VARS, PHASE_DISPLAY, ...
    ALPHA, STRICT_PARAMETRIC_RULE, NcomparisonsHolm);

[Desc_SE, Omni_SE, Pair_SE, LM_SE] = run_raw_value_analysis( ...
    T_SE, 'SampEn', GROUPS, GROUP_DISPLAY, ...
    PANELS, PANEL_DISPLAY, PHASE_VARS, PHASE_DISPLAY, ...
    ALPHA, STRICT_PARAMETRIC_RULE, NcomparisonsHolm);

Desc = [Desc_SWA; Desc_SE];
Omni = [Omni_SWA; Omni_SE];
Pair = [Pair_SWA; Pair_SE];
LMR  = [LM_SWA; LM_SE];

writetable(Desc, fullfile(OUTDIR,'REPORT_DESCRIPTIVES.csv'));
writetable(Omni, fullfile(OUTDIR,'REPORT_OMNIBUS.csv'));
writetable(Pair, fullfile(OUTDIR,'REPORT_PAIRWISE.csv'));
writetable(LMR,  fullfile(OUTDIR,'REPORT_LM.csv'));


%% ===================== PLOTS ============================================
plot_combined_box_raw_values( ...
    T_SWA, T_SE, Pair, Omni, ...
    GROUPS, GROUP_DISPLAY, PANELS, PANEL_DISPLAY, ...
    PHASE_VARS, PHASE_DISPLAY, OUTDIR);

fprintf('\n✅ DONE. Outputs saved in:\n%s\n', OUTDIR);

%% ========================================================================
%% =============================== FUNCTIONS ===============================
%% ========================================================================

function T = load_table_from_saved_mat(matPath)
fprintf('\n[LOAD] %s\n', matPath);

S = load(matPath);

OUT = [];
if isfield(S,'OUT_SWA')
    OUT = S.OUT_SWA;
elseif isfield(S,'OUT_SE')
    OUT = S.OUT_SE;
elseif isfield(S,'OUT')
    OUT = S.OUT;
else
    fn = fieldnames(S);
    for i = 1:numel(fn)
        if isstruct(S.(fn{i})) && isfield(S.(fn{i}),'Table') && istable(S.(fn{i}).Table)
            OUT = S.(fn{i});
            break;
        end
    end
end

if isempty(OUT)
    error('No output struct with Table field found in %s', matPath);
end

if ~isfield(OUT,'Table') || ~istable(OUT.Table)
    error('OUT.Table missing or invalid in %s', matPath);
end

T = OUT.Table;

req = {'Panel','Group','Id','Early','Late'};
for i = 1:numel(req)
    if ~ismember(req{i}, T.Properties.VariableNames)
        error('Missing column "%s" in %s', req{i}, matPath);
    end
end

if isstring(T.Panel), T.Panel = cellstr(T.Panel); end
if isstring(T.Id),    T.Id    = cellstr(T.Id);    end
if iscategorical(T.Group), T.Group = cellstr(T.Group); end

fprintf('  Loaded table: rows=%d | cols=%d\n', height(T), width(T));
disp('  Variables found:');
disp(T.Properties.VariableNames(:)');
end

% ========================================================================
function [RepDesc, RepOmni, RepPair, RepLM] = run_raw_value_analysis( ...
    T, metricName, GROUPS, GROUP_DISPLAY, ...
    PANELS, PANEL_DISPLAY, PHASE_VARS, PHASE_DISPLAY, ...
    ALPHA, STRICT_PARAMETRIC_RULE, NcomparisonsHolm)

RepDesc = table();
RepOmni = table();
RepPair = table();
RepLM   = table();

pairs = nchoosek(1:numel(GROUPS), 2);

fprintf('\n============================================================\n');
fprintf(' ANALYSIS FOR METRIC: %s\n', metricName);
fprintf('============================================================\n');

for ph = 1:numel(PHASE_VARS)
    phaseVar  = PHASE_VARS{ph};
    phaseDisp = PHASE_DISPLAY{ph};

    fprintf('\n============================================================\n');
    fprintf(' PHASE: %s\n', phaseDisp);
    fprintf('============================================================\n');

    for pi = 1:numel(PANELS)
        pan     = PANELS{pi};
        panDisp = PANEL_DISPLAY{pi};

        fprintf('\n------------------------------------------------------------\n');
        fprintf(' PANEL: %s | PHASE: %s\n', panDisp, phaseDisp);
        fprintf('------------------------------------------------------------\n');

        Tp = T(strcmpi(T.Panel, pan), :);
        Tp.Group = categorical(Tp.Group, GROUPS, GROUPS);

        yAll = Tp.(phaseVar);
        ok = isfinite(yAll) & ~isundefined(Tp.Group);
        Tp = Tp(ok,:);

        fprintf('Valid rows after cleanup: %d\n', height(Tp));

        if height(Tp) < 6
            fprintf('SKIP: too few valid rows.\n');
            continue;
        end

        %% ---------- descriptives + normality by group ----------
        groupNormal = nan(numel(GROUPS),1);
        groupN      = nan(numel(GROUPS),1);

        for gi = 1:numel(GROUPS)
            g = GROUPS{gi};
            y = Tp.(phaseVar)(Tp.Group == g);
            y = y(isfinite(y));

            n  = numel(y);
            mu = mean(y,'omitnan');
            sd = std(y,0,'omitnan');
            md = median(y,'omitnan');
            iq = iqr(y);

            [ciLo, ciHi] = ci_mean_t(y, ALPHA);
            [isNormal, pNorm, testUsed] = normality_check_vector(y, ALPHA);

            groupNormal(gi) = isNormal;
            groupN(gi)      = n;

            fprintf('  Group %s | N=%d | mean=%.4g | SD=%.4g | median=%.4g | IQR=%.4g | %s p=%.4g | normal=%s\n', ...
                g, n, mu, sd, md, iq, testUsed, pNorm, logical2str(isNormal));

            RepDesc = [RepDesc; table( ...
                {metricName}, {phaseVar}, {phaseDisp}, {pan}, {panDisp}, {g}, ...
                n, mu, sd, md, iq, ciLo, ciHi, pNorm, isNormal, {testUsed}, ...
                'VariableNames', {'Metric','PhaseVar','PhaseDisplay','Panel','PanelDisplay','Group', ...
                'N','Mean','SD','Median','IQR','CImean_Lo','CImean_Hi','pNormality','NormalityOK','NormalityTest'})]; %#ok<AGROW>
        end

        %% ---------- variance homogeneity ----------
        [pLevene, equalVarOK, leveneTestName] = homogeneity_check(Tp.(phaseVar), Tp.Group, ALPHA);
        fprintf('  Variance homogeneity: %s p=%.4g | equalVarOK=%s\n', ...
            leveneTestName, pLevene, logical2str(equalVarOK));

        %% ---------- residual normality from LM ----------
        lmTmp = fitlm(Tp, sprintf('%s ~ Group', phaseVar));
        resid = lmTmp.Residuals.Raw;
        [residNormal, pResidNorm, residNormTest] = normality_check_vector(resid, ALPHA);

        fprintf('  LM residual normality: %s p=%.4g | normal=%s\n', ...
            residNormTest, pResidNorm, logical2str(residNormal));

        %% ---------- choose omnibus ----------
        allGroupsNormal = all(groupNormal(isfinite(groupNormal)) == 1);
        enoughGroupsForParametric = sum(groupN >= 3) == numel(GROUPS);

        if STRICT_PARAMETRIC_RULE
            useAnova = allGroupsNormal && equalVarOK && residNormal && enoughGroupsForParametric;
        else
            useAnova = allGroupsNormal && enoughGroupsForParametric;
        end

        fprintf('  Decision inputs -> allGroupsNormal=%s | equalVarOK=%s | residNormal=%s | enoughN=%s\n', ...
            logical2str(allGroupsNormal), logical2str(equalVarOK), logical2str(residNormal), logical2str(enoughGroupsForParametric));

        %% ========================== OMNIBUS ===============================
        if useAnova
            fprintf('  ==> Selected OMNIBUS: ONE-WAY ANOVA\n');

            [pOmni, tblAnova, ~] = anova1(Tp.(phaseVar), Tp.Group, 'off');

            SS_between = tblAnova{2,2};
            SS_within  = tblAnova{3,2};
            df1        = tblAnova{2,3};
            df2        = tblAnova{3,3};
            Fval       = tblAnova{2,5};

            eta2   = SS_between / (SS_between + SS_within);
            omega2 = (SS_between - df1*(SS_within/df2)) / (SS_between + SS_within + (SS_within/df2));
            omega2 = max(0, omega2);

            fprintf('     ANOVA: F(%g,%g)=%.4f | p=%.4g | eta2=%.4f | omega2=%.4f\n', ...
                df1, df2, Fval, pOmni, eta2, omega2);

            omnibusType = 'ANOVA';

            RepOmni = [RepOmni; table( ...
                {metricName}, {phaseVar}, {phaseDisp}, {pan}, {panDisp}, height(Tp), ...
                allGroupsNormal, equalVarOK, residNormal, pResidNorm, pLevene, ...
                {omnibusType}, pOmni, Fval, df1, df2, eta2, omega2, NaN, NcomparisonsHolm, ...
                'VariableNames', {'Metric','PhaseVar','PhaseDisplay','Panel','PanelDisplay','Ntotal', ...
                'AllGroupsNormal','EqualVarianceOK','ResidualNormalOK','pResidualNormality','pVarianceHomogeneity', ...
                'OmnibusTest','pOmnibus','Stat','DF1','DF2','Eta2','Omega2','Epsilon2','NcomparisonsHolm'})]; %#ok<AGROW>

            %% ---------- pairwise after ANOVA ----------
            pRaw = nan(size(pairs,1),1);
            tmp  = cell(size(pairs,1),1);

            for k = 1:size(pairs,1)
                g1 = GROUPS{pairs(k,1)};
                g2 = GROUPS{pairs(k,2)};

                y1 = Tp.(phaseVar)(Tp.Group == g1); y1 = y1(isfinite(y1));
                y2 = Tp.(phaseVar)(Tp.Group == g2); y2 = y2(isfinite(y2));

                mean1 = mean(y1,'omitnan');
                mean2 = mean(y2,'omitnan');
                sd1   = std(y1,0,'omitnan');
                sd2   = std(y2,0,'omitnan');

                [~,p,ci,stats] = ttest2(y1, y2, 'Vartype', tern(equalVarOK,'equal','unequal'));

                mdiff   = mean2 - mean1;
                hedgesg = hedges_g_independent(y2, y1);

                pRaw(k) = p;
                tmp{k} = {metricName,phaseVar,phaseDisp,pan,panDisp,omnibusType,g1,g2, ...
                    numel(y1),numel(y2), ...
                    mean1,mean2,sd1,sd2, ...
                    mdiff,ci(1),ci(2),stats.tstat,stats.df,p,NaN,hedgesg,'ttest2',NcomparisonsHolm};

                fprintf('     Pair %s vs %s | t=%.4f | df=%.3f | p=%.4g | meanDiff=%.4g | CI95%%=[%.4g, %.4g] | Hedges g=%.4f\n', ...
                    g1, g2, stats.tstat, stats.df, p, mdiff, ci(1), ci(2), hedgesg);
            end

            pAdj = apply_holm_vector_n(pRaw, NcomparisonsHolm);

            for k = 1:size(pairs,1)
                tmp{k}{21} = pAdj(k);
                RepPair = [RepPair; cell2table(tmp{k}, 'VariableNames', ...
                    {'Metric','PhaseVar','PhaseDisplay','Panel','PanelDisplay','OmnibusType','Group1','Group2', ...
                     'N1','N2','Mean1','Mean2','SD1','SD2', ...
                     'EffectEstimate','CI_Lo','CI_Hi','Stat','DF','pRaw','pHolm','EffectSize','PairwiseTest','NcomparisonsHolm'})]; %#ok<AGROW>
            end

        else
            fprintf('  ==> Selected OMNIBUS: KRUSKAL-WALLIS\n');

            [pOmni, tblKW, ~] = kruskalwallis(Tp.(phaseVar), Tp.Group, 'off');

            chi2 = tblKW{2,5};
            df1  = tblKW{2,3};
            nTot = height(Tp);

            epsilon2 = (chi2 - df1 + 1) / (nTot - df1);
            epsilon2 = max(0, epsilon2);

            fprintf('     KW: chi2(%g)=%.4f | p=%.4g | epsilon2=%.4f\n', ...
                df1, chi2, pOmni, epsilon2);

            omnibusType = 'KruskalWallis';

            RepOmni = [RepOmni; table( ...
                {metricName}, {phaseVar}, {phaseDisp}, {pan}, {panDisp}, height(Tp), ...
                allGroupsNormal, equalVarOK, residNormal, pResidNorm, pLevene, ...
                {omnibusType}, pOmni, chi2, df1, NaN, NaN, NaN, epsilon2, NcomparisonsHolm, ...
                'VariableNames', {'Metric','PhaseVar','PhaseDisplay','Panel','PanelDisplay','Ntotal', ...
                'AllGroupsNormal','EqualVarianceOK','ResidualNormalOK','pResidualNormality','pVarianceHomogeneity', ...
                'OmnibusTest','pOmnibus','Stat','DF1','DF2','Eta2','Omega2','Epsilon2','NcomparisonsHolm'})]; %#ok<AGROW>

            %% ---------- pairwise after KW ----------
            pRaw = nan(size(pairs,1),1);
            tmp  = cell(size(pairs,1),1);

            for k = 1:size(pairs,1)
                g1 = GROUPS{pairs(k,1)};
                g2 = GROUPS{pairs(k,2)};

                y1 = Tp.(phaseVar)(Tp.Group == g1); y1 = y1(isfinite(y1));
                y2 = Tp.(phaseVar)(Tp.Group == g2); y2 = y2(isfinite(y2));

                p = ranksum(y1, y2);

                mean1 = mean(y1,'omitnan');
                mean2 = mean(y2,'omitnan');
                sd1   = std(y1,0,'omitnan');
                sd2   = std(y2,0,'omitnan');

                [hl, ciLo, ciHi] = hodges_lehmann_ci(y2, y1, ALPHA, 4000);
                rb = rank_biserial_from_groups(y2, y1);

                pRaw(k) = p;
                tmp{k} = {metricName,phaseVar,phaseDisp,pan,panDisp,omnibusType,g1,g2, ...
                    numel(y1),numel(y2), ...
                    mean1,mean2,sd1,sd2, ...
                    hl,ciLo,ciHi,NaN,NaN,p,NaN,rb,'ranksum',NcomparisonsHolm};

                fprintf('     Pair %s vs %s | ranksum p=%.4g | HL shift=%.4g | CI95%%=[%.4g, %.4g] | rank-biserial=%.4f\n', ...
                    g1, g2, p, hl, ciLo, ciHi, rb);
            end

            pAdj = apply_holm_vector_n(pRaw, NcomparisonsHolm);

            for k = 1:size(pairs,1)
                tmp{k}{21} = pAdj(k);
                RepPair = [RepPair; cell2table(tmp{k}, 'VariableNames', ...
                    {'Metric','PhaseVar','PhaseDisplay','Panel','PanelDisplay','OmnibusType','Group1','Group2', ...
                     'N1','N2','Mean1','Mean2','SD1','SD2', ...
                     'EffectEstimate','CI_Lo','CI_Hi','Stat','DF','pRaw','pHolm','EffectSize','PairwiseTest','NcomparisonsHolm'})]; %#ok<AGROW>
            end
        end

        %% =========================== LM CHECK =============================
        lm = fitlm(Tp, sprintf('%s ~ Group', phaseVar));
        a  = anova(lm, 'summary');

        idxGroup = find(strcmp(a.Properties.RowNames,'Group'),1,'first');
        if isempty(idxGroup)
            pGroup = NaN; FGroup = NaN; dfGroup = NaN;
        else
            pGroup  = a.pValue(idxGroup);
            FGroup  = a.F(idxGroup);
            dfGroup = a.DF(idxGroup);
        end

        fprintf('  LM CHECK: %s ~ Group\n', phaseVar);
        fprintf('     Group effect: F=%.4f | df=%g | p=%.4g | MSE=%.4g | R2=%.4f | adjR2=%.4f\n', ...
            FGroup, dfGroup, pGroup, lm.MSE, lm.Rsquared.Ordinary, lm.Rsquared.Adjusted);

        coefTbl = lm.Coefficients;
        for rr = 1:height(coefTbl)
            fprintf('     Coef %-18s | beta=% .5g | SE=%.5g | t=% .4f | p=%.4g\n', ...
                coefTbl.Properties.RowNames{rr}, ...
                coefTbl.Estimate(rr), coefTbl.SE(rr), coefTbl.tStat(rr), coefTbl.pValue(rr));
        end

        RepLM = [RepLM; table( ...
            {metricName}, {phaseVar}, {phaseDisp}, {pan}, {panDisp}, height(Tp), ...
            pGroup, FGroup, dfGroup, lm.DFE, lm.MSE, lm.Rsquared.Ordinary, lm.Rsquared.Adjusted, pResidNorm, residNormal, ...
            'VariableNames', {'Metric','PhaseVar','PhaseDisplay','Panel','PanelDisplay','Ntotal', ...
            'pGroup_LM','F_Group_LM','DF_Group_LM','DFE_LM','MSE_LM','R2','AdjR2','pResidualNormality','ResidualNormalOK'})]; %#ok<AGROW>
    end
end
end

% ========================================================================
function [isNormal, pVal, testName] = normality_check_vector(x, alpha)
x = x(:);
x = x(isfinite(x));

if numel(x) < 3
    isNormal = NaN;
    pVal = NaN;
    testName = 'insufficientN';
    return;
end

if exist('swtest','file') == 2
    try
        [h,pVal] = swtest(x, alpha);
        isNormal = double(h == 0);
        testName = 'swtest';
        return;
    catch
    end
end

try
    [h,pVal] = lillietest(x, 'Alpha', alpha);
    isNormal = double(h == 0);
    testName = 'lillietest';
catch
    isNormal = NaN;
    pVal = NaN;
    testName = 'normalityTestFailed';
end
end

% ========================================================================
function [pVal, ok, testName] = homogeneity_check(y, g, alpha)
y = y(:);
okRows = isfinite(y) & ~isundefined(g);
y = y(okRows);
g = g(okRows);

try
    pVal = vartestn(y, g, 'TestType','LeveneAbsolute', 'Display','off');
    ok = double(pVal >= alpha);
    testName = 'LeveneAbsolute';
catch
    pVal = NaN;
    ok = NaN;
    testName = 'LeveneFailed';
end
end

% ========================================================================
function [ciLo, ciHi] = ci_mean_t(x, alpha)
x = x(:);
x = x(isfinite(x));
n = numel(x);

if n < 2
    ciLo = NaN;
    ciHi = NaN;
    return;
end

m = mean(x);
s = std(x,0);
se = s / sqrt(n);
tcrit = tinv(1 - alpha/2, n-1);

ciLo = m - tcrit*se;
ciHi = m + tcrit*se;
end

% ========================================================================
function pAdj = apply_holm_vector_n(p, NcomparisonsHolm)
pAdj = nan(size(p));
ok = isfinite(p);
pv = p(ok);

mObs = numel(pv);
if mObs == 0
    return;
end

m = max(NcomparisonsHolm, mObs);

[ps, ord] = sort(pv, 'ascend');
adj = nan(size(ps));

for i = 1:mObs
    adj(i) = (m - i + 1) * ps(i);
end
for i = mObs-1:-1:1
    adj(i) = min(adj(i), adj(i+1));
end
adj = min(adj, 1);

back = nan(size(pv));
back(ord) = adj;
pAdj(ok) = back;
end

% ========================================================================
function g = hedges_g_independent(x1, x2)
x1 = x1(:); x2 = x2(:);
x1 = x1(isfinite(x1));
x2 = x2(isfinite(x2));

n1 = numel(x1);
n2 = numel(x2);

if n1 < 2 || n2 < 2
    g = NaN;
    return;
end

m1 = mean(x1);
m2 = mean(x2);
s1 = std(x1,0);
s2 = std(x2,0);

sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2));
if sp <= 0 || ~isfinite(sp)
    g = NaN;
    return;
end

d = (m1 - m2) / sp;
J = 1 - 3/(4*(n1+n2)-9);
g = d * J;
end

% ========================================================================
function rb = rank_biserial_from_groups(x1, x2)
x1 = x1(:); x2 = x2(:);
x1 = x1(isfinite(x1));
x2 = x2(isfinite(x2));

n1 = numel(x1);
n2 = numel(x2);

if n1 == 0 || n2 == 0
    rb = NaN;
    return;
end

allVals = [x1; x2];
ranks = tiedrank(allVals);
R1 = sum(ranks(1:n1));
U1 = R1 - n1*(n1+1)/2;

rb = 2*U1/(n1*n2) - 1;
end

% ========================================================================
function [hl, ciLo, ciHi] = hodges_lehmann_ci(x1, x2, alpha, nBoot)
if nargin < 3 || isempty(alpha), alpha = 0.05; end
if nargin < 4 || isempty(nBoot), nBoot = 4000; end

x1 = x1(:);
x2 = x2(:);
x1 = x1(isfinite(x1));
x2 = x2(isfinite(x2));

if isempty(x1) || isempty(x2)
    hl   = NaN;
    ciLo = NaN;
    ciHi = NaN;
    return;
end

D  = x1 - x2';
hl = median(D(:), 'omitnan');

n1 = numel(x1);
n2 = numel(x2);
bootStat = nan(nBoot,1);

for b = 1:nBoot
    idx1 = randi(n1, [n1 1]);
    idx2 = randi(n2, [n2 1]);
    xb1 = x1(idx1);
    xb2 = x2(idx2);
    Db = xb1 - xb2';
    bootStat(b) = median(Db(:), 'omitnan');
end

ciLo = prctile(bootStat, 100*(alpha/2));
ciHi = prctile(bootStat, 100*(1-alpha/2));
end

% ========================================================================
function out = tern(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

% ========================================================================
function s = logical2str(x)
if isnan(x)
    s = 'NaN';
elseif x
    s = 'YES';
else
    s = 'NO';
end
end

% ========================================================================
function plot_combined_box_raw_values(T_SWA, T_SE, Pair, Omni, GROUPS, GROUP_DISPLAY, PANELS, PANEL_DISPLAY, PHASE_VARS, PHASE_DISPLAY, OUTDIR)

C = [0.70 0.82 0.95;
     0.98 0.82 0.67;
     0.76 0.90 0.76];

metrics = {'SWA_classic','SampEn'};
Tcell   = {T_SWA, T_SE};

nRows = numel(metrics) * numel(PHASE_VARS);
nCols = numel(PANELS);

% y-limits per row = metric × phase
ylims = nan(nRows,2);
rowCounter = 0;

for m = 1:numel(metrics)
    Tm = Tcell{m};
    for ph = 1:numel(PHASE_VARS)
        rowCounter = rowCounter + 1;
        phaseVar = PHASE_VARS{ph};

        yAll = [];
        for pi = 1:numel(PANELS)
            Tp = Tm(strcmpi(Tm.Panel, PANELS{pi}), :);
            y = Tp.(phaseVar);
            y = y(isfinite(y));
            yAll = [yAll; y(:)]; %#ok<AGROW>
        end

        if isempty(yAll)
            ylims(rowCounter,:) = [-1 1];
        else
            lo = prctile(yAll,2);
            hi = prctile(yAll,98);
            pad = 0.08*(hi-lo + eps);
            ylims(rowCounter,:) = [lo-pad, hi+pad];
        end
    end
end

fig = figure('Color','w','Position',[80 80 1550 1100]);
tlo = tiledlayout(fig, nRows, nCols, 'Padding','compact', 'TileSpacing','compact');

rowCounter = 0;
for m = 1:numel(metrics)
    metricName = metrics{m};
    Tm = Tcell{m};

    for ph = 1:numel(PHASE_VARS)
        rowCounter = rowCounter + 1;

        phaseVar  = PHASE_VARS{ph};
        phaseDisp = PHASE_DISPLAY{ph};

        for pi = 1:numel(PANELS)
            pan = PANELS{pi};
            ax = nexttile(tlo, (rowCounter-1)*nCols + pi);
            hold(ax,'on');

            Tp = Tm(strcmpi(Tm.Panel, pan), :);
            Tp.Group = categorical(Tp.Group, GROUPS, GROUPS);

            y = Tp.(phaseVar);
            g = Tp.Group;
            ok = isfinite(y) & ~isundefined(g);
            y = y(ok);
            g = g(ok);

            set(ax,'Color','w');
            grid(ax,'off');
            box(ax,'off');
            ax.LineWidth = 1.0;
            ax.FontSize = 10;

            for gi = 1:numel(GROUPS)
                yy = y(g==GROUPS{gi});
                if isempty(yy), continue; end
                xx = gi*ones(size(yy));
                bc = boxchart(ax, xx, yy);
                bc.BoxWidth = 0.55;
                bc.MarkerStyle = 'none';
                bc.BoxFaceColor = C(gi,:);
                bc.BoxFaceAlpha = 0.45;
                bc.LineWidth = 1.2;
            end

            xNum = double(g);
            jitter = (rand(size(xNum))-0.5)*0.22;
            scatter(ax, xNum + jitter, y, 16, 'filled', ...
                'MarkerFaceAlpha',0.30, 'MarkerEdgeAlpha',0);

            ax.XTick = 1:numel(GROUPS);
            ax.XTickLabel = GROUP_DISPLAY;
            ax.TickLabelInterpreter = 'tex';
            ylim(ax, ylims(rowCounter,:));

            if rowCounter == 1
                title(ax, PANEL_DISPLAY{pi}, 'FontWeight','bold');
            end

            if pi == 1
                ylabel(ax, sprintf('%s (%s)', metricName, phaseDisp), 'FontWeight','bold');
            end

            idx = strcmp(Omni.Metric, metricName) & strcmp(Omni.Panel, pan) & strcmp(Omni.PhaseVar, phaseVar);
            if any(idx)
                row = find(idx,1,'first');
                pOmni = Omni.pOmnibus(row);
                testN = Omni.OmnibusTest{row};

                if strcmpi(testN,'ANOVA')
                    txt = sprintf('AN p=%.3g | \\eta^2=%.2f', pOmni, Omni.Eta2(row));
                else
                    txt = sprintf('KW p=%.3g | \\epsilon^2=%.2f', pOmni, Omni.Epsilon2(row));
                end

                text(ax, 0.03, 0.97, txt, ...
                    'Units','normalized', ...
                    'HorizontalAlignment','left', ...
                    'VerticalAlignment','top', ...
                    'BackgroundColor','w', ...
                    'Margin',3);
            end

            add_sig_segments_group_phase(ax, Pair, metricName, phaseVar, pan, GROUPS);
        end
    end
end

base = fullfile(OUTDIR, 'FIG_COMBINED_RAW_EARLY_LATE_GROUPONLY');
exportgraphics(fig,[base '.png'],'Resolution',300);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');

fprintf('Saved combined figure:\n%s\n', base);
end

% ========================================================================
function add_sig_segments_group_phase(ax, Pair, metricName, phaseVar, panelKey, GROUPS)

if isempty(Pair) || ~istable(Pair)
    return;
end

idx = strcmp(Pair.Metric, metricName) & ...
      strcmp(Pair.PhaseVar, phaseVar) & ...
      strcmp(Pair.Panel, panelKey) & ...
      isfinite(Pair.pHolm);

P = Pair(idx,:);
if isempty(P), return; end

P = P(P.pHolm < 0.05, :);
if isempty(P), return; end

yl   = ylim(ax);
span = yl(2) - yl(1);

step  = 0.06 * span;
hTick = 0.02 * span;
yBase = yl(2);
kStack = 0;

[~,ord] = sort(P.pHolm, 'ascend');
P = P(ord,:);

for i = 1:height(P)
    g1 = P.Group1{i};
    g2 = P.Group2{i};

    x1 = find(strcmp(GROUPS, g1));
    x2 = find(strcmp(GROUPS, g2));
    if isempty(x1) || isempty(x2)
        continue;
    end

    kStack = kStack + 1;
    y = yBase + kStack * step;

    if kStack == 1
        y = y - 0.025 * span;
    end

    plot(ax, [x1 x1 x2 x2], [y-hTick y y y-hTick], 'k-', 'LineWidth', 1.2);

    p = P.pHolm(i);
    if p < 1e-3
        stars = '***';
    elseif p < 1e-2
        stars = '**';
    else
        stars = '*';
    end

    if kStack == 1
        starOffset = 0.004 * span;
    else
        starOffset = 0.010 * span;
    end

    text(ax, mean([x1 x2]), y + starOffset, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold', ...
        'FontSize', 11);
end

ylim(ax, [yl(1), yBase + (kStack + 2) * step]);
end