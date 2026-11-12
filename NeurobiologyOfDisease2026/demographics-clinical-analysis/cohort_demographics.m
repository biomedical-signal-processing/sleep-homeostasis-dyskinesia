%% ========================================================================
% Table 1 – Demographics/Clinical re-analysis (CTL vs ADV vs DYS)
% IMPORTANT: some variables are PD-only -> test ADV vs DYS only
%
% INPUT:
%   db_path    : db_annals_theta_CLEAN.xlsx
%   swtest_dir : folder containing swtest.m (Shapiro–Wilk)
%
% OUTPUT:
%   Table1_Recap.xlsx with sheets:
%       Table1_formatted, Assumptions, Omnibus, Run_log
%
% Tests:
%   Normality: Shapiro–Wilk (swtest) within each included group
%   Homoscedasticity:
%     - 3-group: Levene (vartestn, LeveneAbsolute)
%     - 2-group: F-test on variances (vartest2)
%   Omnibus:
%     - 3-group:
%         ANOVA if all normal & Levene p>0.05
%         Welch ANOVA if all normal & Levene p<=0.05
%         Kruskal–Wallis otherwise
%     - 2-group (ADV vs DYS):
%         t-test if both normal & equal variances
%         Welch t-test if both normal & unequal variances
%         ranksum (Mann–Whitney) otherwise
%   Sex:
%     - 3-group: Chi-square (crosstab), effect size Cramer's V
%
% Effect sizes:
%   - ANOVA: eta^2
%   - Welch ANOVA: eta^2 (approx)
%   - Kruskal–Wallis: epsilon^2
%   - t-test/Welch t-test: Cohen's d (pooled SD; for Welch still reported as d)
%   - ranksum: rank-biserial r
%   - Chi-square: Cramer's V
%
% NOTE: ESS, PSQI are available in CTL -> included in 3-group tests.
% ========================================================================

clear; clc; close all;

%% ---------------- PATHS (EDIT) ----------------
db_path    = 'C:\Users\EOC\OneDrive\script matlab\Annals_entropy\database\db_annals_theta_CLEAN.xlsx';
swtest_dir = 'C:\Users\EOC\OneDrive\script matlab\appoggio\swtest';
addpath(swtest_dir);

out_xlsx = fullfile(fileparts(db_path), 'Table1_Recap.xlsx');

%% ---------------- LOAD DB ----------------
opts = detectImportOptions(db_path);
opts.VariableNamingRule = 'preserve';
T = readtable(db_path, opts);

% Normalize Group
T.Group = upper(strtrim(string(T.Group)));
keepGroups = ["CTL","ADV","DYS"];
T = T(ismember(T.Group, keepGroups), :);

% Normalize Sex to 'F'/'M'
T.Sex = upper(strtrim(string(T.Sex)));

%% ---------------- VARIABLES (Table 1) ----------------
% name in DB | label | type | scope
% scope: 'ALL' (CTL/ADV/DYS) or 'PD' (ADV/DYS only)
vars = {
    'Age',             'Age, yr',              'numeric',     'ALL';
    'Sex',             'Gender, % (F/M)',      'categorical', 'ALL';

    'DiseaseDuration', 'Disease duration, yr', 'numeric',     'PD';
    'H&Y',             'H&Y',                  'numeric',     'PD';
    'UPDRS_I',         'MDS-UPDRS I',          'numeric',     'PD';
    'UPDRS_II',        'MDS-UPDRS II',         'numeric',     'PD';
    'UPDRS_III',       'MDS-UPDRS III',        'numeric',     'PD';
    'UPDRS_IV',        'MDS-UPDRS IV',         'numeric',     'PD';
    'UPDRS_Tot',       'MDS-UPDRS Tot',        'numeric',     'PD';
    'AIMS',            'AIMS',                 'numeric',     'PD';
    'LEDD',            'LEDD, mg',             'numeric',     'PD';

    'ESS',             'ESS',                  'numeric',     'ALL';
    'PSQI',            'PSQI',                 'numeric',     'ALL';
};

%% ---------------- PREP OUTPUTS ----------------
nV = size(vars,1);
outVar   = strings(nV,1);
col_CTL  = strings(nV,1);
col_ADV  = strings(nV,1);
col_DYS  = strings(nV,1);

testName = strings(nV,1);
statName = strings(nV,1);
df1      = nan(nV,1);
df2      = nan(nV,1);
pval     = nan(nV,1);
effName  = strings(nV,1);
effVal   = nan(nV,1);

assump_rows = {}; % Variable, Group, N, Shapiro_p, Normality_pass, Levene/Ftest_p, Homosced_pass, Notes
omni_rows   = {}; % Variable, Scope, Test, StatisticName, Statistic, df1, df2, p, EffectSizeName, EffectSizeValue, Notes

G = T.Group;
idx_CTL = (G=="CTL");
idx_ADV = (G=="ADV");
idx_DYS = (G=="DYS");

n_CTL = sum(idx_CTL);
n_ADV = sum(idx_ADV);
n_DYS = sum(idx_DYS);

%% ---------------- LOOP ----------------
for i = 1:nV

    vname  = vars{i,1};
    vlabel = vars{i,2};
    vtype  = vars{i,3};
    scope  = vars{i,4};

    outVar(i) = vlabel;

    % included groups
    if scope == "ALL"
        use_CTL = true; use_ADV = true; use_DYS = true;
    else
        use_CTL = false; use_ADV = true; use_DYS = true;
    end

    % ---------------- CATEGORICAL: Sex (ALL only) ----------------
    if strcmp(vtype,'categorical')
        x = upper(strtrim(string(T.(vname))));

        % counts for display
        cCTL = x(idx_CTL); cADV = x(idx_ADV); cDYS = x(idx_DYS);

        nF_CTL = sum(cCTL=="F"); nM_CTL = sum(cCTL=="M");
        nF_ADV = sum(cADV=="F"); nM_ADV = sum(cADV=="M");
        nF_DYS = sum(cDYS=="F"); nM_DYS = sum(cDYS=="M");

        col_CTL(i) = sprintf('%d/%d (%.1f/%.1f)', nF_CTL, nM_CTL, 100*nF_CTL/max(n_CTL,1), 100*nM_CTL/max(n_CTL,1));
        col_ADV(i) = sprintf('%d/%d (%.1f/%.1f)', nF_ADV, nM_ADV, 100*nF_ADV/max(n_ADV,1), 100*nM_ADV/max(n_ADV,1));
        col_DYS(i) = sprintf('%d/%d (%.1f/%.1f)', nF_DYS, nM_DYS, 100*nF_DYS/max(n_DYS,1), 100*nM_DYS/max(n_DYS,1));

        [tbl,~,~,~] = crosstab(G, x);

        expected = sum(tbl,2)*sum(tbl,1)/sum(tbl(:));
        chi2 = sum((tbl-expected).^2 ./ max(expected,eps), 'all');
        df_chi = (size(tbl,1)-1)*(size(tbl,2)-1);
        p = 1-chi2cdf(chi2, df_chi);

        testName(i) = "Chi-square";
        statName(i) = "chi2";
        df1(i) = df_chi;
        df2(i) = NaN;
        pval(i)= p;

        nTot = sum(tbl(:));
        kmin = min(size(tbl,1)-1, size(tbl,2)-1);
        V = sqrt(chi2 / (nTot * max(kmin,1)));
        effName(i) = "Cramer's V";
        effVal(i)  = V;

        omni_rows(end+1,:) = {vlabel, "ALL", "Chi-square", "chi2", chi2, df_chi, NaN, p, "Cramer's V", V, ""}; %#ok<SAGROW>
        continue;
    end

    % ---------------- NUMERIC ----------------
    y_all = double(T.(vname));
    g_all = string(G);

    % Filter by scope
    if scope == "ALL"
        keep = ~isnan(y_all) & (g_all=="CTL" | g_all=="ADV" | g_all=="DYS");
    else
        keep = ~isnan(y_all) & (g_all=="ADV" | g_all=="DYS");
    end
    y = y_all(keep);
    g = g_all(keep);

    y_CTL = y(g=="CTL");
    y_ADV = y(g=="ADV");
    y_DYS = y(g=="DYS");

    % Descriptives for display
    if use_CTL
        mu_CTL = mean(y_CTL,'omitnan'); sd_CTL = std(y_CTL,0,'omitnan');
        col_CTL(i) = sprintf('%.2f ± %.2f', mu_CTL, sd_CTL);
    else
        col_CTL(i) = "—";
    end

    mu_ADV = mean(y_ADV,'omitnan'); sd_ADV = std(y_ADV,0,'omitnan');
    mu_DYS = mean(y_DYS,'omitnan'); sd_DYS = std(y_DYS,0,'omitnan');
    col_ADV(i) = sprintf('%.2f ± %.2f', mu_ADV, sd_ADV);
    col_DYS(i) = sprintf('%.2f ± %.2f', mu_DYS, sd_DYS);

    % Shapiro per gruppo (swtest)
    pSW_CTL = NaN; pSW_ADV = NaN; pSW_DYS = NaN;
    try
        if use_CTL && numel(y_CTL) >= 3, [~, pSW_CTL] = swtest(y_CTL); end
        if numel(y_ADV) >= 3, [~, pSW_ADV] = swtest(y_ADV); end
        if numel(y_DYS) >= 3, [~, pSW_DYS] = swtest(y_DYS); end
    catch
    end

    % Homoscedasticity test
    pHom = NaN;
    homo = NaN;

    if scope == "ALL"
        % Levene across groups
        try
            pHom = vartestn(y, g, 'TestType','LeveneAbsolute', 'Display','off');
            homo = (pHom > 0.05);
        catch
        end
        % Normality criterion: all included groups normal
        allNormal = true;
        if use_CTL, allNormal = allNormal & (isnan(pSW_CTL) | pSW_CTL>0.05); end
        allNormal = allNormal & (isnan(pSW_ADV) | pSW_ADV>0.05);
        allNormal = allNormal & (isnan(pSW_DYS) | pSW_DYS>0.05);

        % Store assumptions rows
        if use_CTL, assump_rows(end+1,:) = {vlabel,"CTL",numel(y_CTL),pSW_CTL,allNormal,pHom,homo,"Levene across groups"}; end %#ok<SAGROW>
        assump_rows(end+1,:) = {vlabel,"ADV",numel(y_ADV),pSW_ADV,allNormal,pHom,homo,"Levene across groups"}; %#ok<SAGROW>
        assump_rows(end+1,:) = {vlabel,"DYS",numel(y_DYS),pSW_DYS,allNormal,pHom,homo,"Levene across groups"}; %#ok<SAGROW>

        % Omnibus choice (3 groups)
        if allNormal && homo
            [p, tbl] = anova1(y, g, 'off');
            F = tbl{2,5};
            df_between = tbl{2,3};
            df_within  = tbl{3,3};

            testName(i) = "ANOVA";
            statName(i) = "F";
            df1(i) = df_between;
            df2(i) = df_within;
            pval(i)= p;

            SSb = tbl{2,2}; SSt = tbl{4,2};
            eta2 = SSb / SSt;

            effName(i) = "eta^2";
            effVal(i)  = eta2;

            omni_rows(end+1,:) = {vlabel,"ALL","ANOVA","F",F,df_between,df_within,p,"eta^2",eta2,""}; %#ok<SAGROW>

        elseif allNormal && ~homo
            [F, dfn, dfd, p] = welch_anova_oneway(y, g);

            testName(i) = "Welch ANOVA";
            statName(i) = "F";
            df1(i) = dfn;
            df2(i) = dfd;
            pval(i)= p;

            % approx eta^2
            grpNames = unique(g);
            grand = mean(y,'omitnan');
            SSb = 0;
            for gg = 1:numel(grpNames)
                yy = y(g==grpNames(gg));
                SSb = SSb + numel(yy) * (mean(yy,'omitnan') - grand).^2;
            end
            SSt = sum((y - grand).^2);
            eta2 = SSb / SSt;

            effName(i) = "eta^2 (approx)";
            effVal(i)  = eta2;

            omni_rows(end+1,:) = {vlabel,"ALL","Welch ANOVA","F",F,dfn,dfd,p,"eta^2 (approx)",eta2,""}; %#ok<SAGROW>

        else
            [p, tbl] = kruskalwallis(y, g, 'off');
            H = tbl{2,5};
            df_kw = tbl{2,3};

            testName(i) = "Kruskal-Wallis";
            statName(i) = "H";
            df1(i) = df_kw;
            df2(i) = NaN;
            pval(i)= p;

            kgrp = numel(unique(g));
            nTot = numel(y);
            eps2 = max(0, (H - kgrp + 1) / (nTot - kgrp));

            effName(i) = "epsilon^2";
            effVal(i)  = eps2;

            omni_rows(end+1,:) = {vlabel,"ALL","Kruskal-Wallis","H",H,df_kw,NaN,p,"epsilon^2",eps2,""}; %#ok<SAGROW>
        end

    else
        % scope == PD (ADV vs DYS only)
        % F-test for equal variances (only if both normal)
        bothNormal = ( (isnan(pSW_ADV) | pSW_ADV>0.05) & (isnan(pSW_DYS) | pSW_DYS>0.05) );

        pHom = NaN; homo = NaN;
        if bothNormal && numel(y_ADV)>=2 && numel(y_DYS)>=2
            try
                [hVar, pHom] = vartest2(y_ADV, y_DYS); %#ok<ASGLU>
                homo = (pHom > 0.05);
            catch
            end
        end

        assump_rows(end+1,:) = {vlabel,"ADV",numel(y_ADV),pSW_ADV,bothNormal,pHom,homo,"F-test var (2 groups)"}; %#ok<SAGROW>
        assump_rows(end+1,:) = {vlabel,"DYS",numel(y_DYS),pSW_DYS,bothNormal,pHom,homo,"F-test var (2 groups)"}; %#ok<SAGROW>

        if bothNormal
            % t-test (equal var) or Welch t-test (unequal var)
            vtype_t = 'equal';
            if ~isnan(homo) && ~homo
                vtype_t = 'unequal';
            end
            [h, p, ci, stats] = ttest2(y_ADV, y_DYS, 'Vartype', vtype_t); %#ok<ASGLU>
            tstat = stats.tstat;
            df_t  = stats.df;

            testName(i) = "t-test (ADV vs DYS)";
            if strcmp(vtype_t,'unequal')
                testName(i) = "Welch t-test (ADV vs DYS)";
            end
            statName(i) = "t";
            df1(i) = df_t; df2(i)=NaN;
            pval(i)= p;

            % Cohen's d
            n1 = numel(y_ADV); n2 = numel(y_DYS);
            s1 = std(y_ADV,0,'omitnan'); s2 = std(y_DYS,0,'omitnan');
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(n1+n2-2,1));
            d = (mean(y_ADV,'omitnan') - mean(y_DYS,'omitnan')) / max(sp,eps);

            effName(i) = "Cohen's d";
            effVal(i)  = d;

            omni_rows(end+1,:) = {vlabel,"PD","ttest2","t",tstat,df_t,NaN,p,"Cohen's d",d,""}; %#ok<SAGROW>

        else
            % Mann-Whitney (ranksum)
            [p, h, stats] = ranksum(y_ADV, y_DYS); %#ok<ASGLU>
            U = stats.ranksum; % MATLAB reports ranksum; U derivation not always direct
            testName(i) = "Mann-Whitney (ADV vs DYS)";
            statName(i) = "U";
            df1(i) = NaN; df2(i)=NaN;
            pval(i)= p;

            % Rank-biserial r from U:
            n1 = numel(y_ADV); n2 = numel(y_DYS);
            U1 = U - n1*(n1+1)/2;
            rrb = 1 - (2*U1)/(n1*n2); % rank-biserial

            effName(i) = "rank-biserial r";
            effVal(i)  = rrb;

            omni_rows(end+1,:) = {vlabel,"PD","ranksum","U",U,NaN,NaN,p,"rank-biserial r",rrb,""}; %#ok<SAGROW>
        end
    end

end

%% ---------------- BUILD TABLES ----------------
Table1_formatted = table(outVar, col_CTL, col_ADV, col_DYS, testName, statName, df1, df2, pval, effName, effVal, ...
    'VariableNames', {'Variable','CTL','ADV','DYS','Test','Stat','df1','df2','p','EffectSize','EffectValue'});

Assumptions = cell2table(assump_rows, 'VariableNames', ...
    {'Variable','Group','N','Shapiro_p','Normality_pass','Homosced_p','Homosced_pass','Notes'});

Omnibus = cell2table(omni_rows, 'VariableNames', ...
    {'Variable','Scope','Test','StatisticName','Statistic','df1','df2','p','EffectSizeName','EffectSizeValue','Notes'});

Run_log = table(string(datetime('now')), string(db_path), string(swtest_dir), ...
    "CTL, ADV, DYS", string(strjoin(string(vars(:,2))', ', ')), ...
    "Generated by Table1_demographics_analysis_PDaware.m", ...
    'VariableNames', {'Timestamp','DB_path','swtest_path','Groups_included','Variables_analysed','Notes'});

%% ---------------- WRITE EXCEL ----------------
if exist(out_xlsx,'file'), delete(out_xlsx); end
writetable(Table1_formatted, out_xlsx, 'Sheet','Table1_formatted');
writetable(Assumptions,      out_xlsx, 'Sheet','Assumptions');
writetable(Omnibus,          out_xlsx, 'Sheet','Omnibus');
writetable(Run_log,          out_xlsx, 'Sheet','Run_log');

fprintf('\nSaved: %s\n', out_xlsx);
