
%% ========================================================================
% SCRIPT: Partial correlations between log-transformed theta power
%         and clinical variables in Parkinson's disease
%
% DESCRIPTION:
% This script performs partial correlation analyses between
% log-transformed theta EEG power (4–8 Hz) and clinical variables
% (LEDD and UDysRS Part IV) in Parkinson's disease patients.
%
% Analysis 1: Theta morning vs LEDD
%   - Groups: EPD, ADV, DYS
%
% Analysis 2: Diurnal theta upscaling vs LEDD
%   - Evening minus morning log-transformed theta power
%   - Groups: EPD, ADV, DYS
%
% Analysis 3: Theta morning vs UDysRS Part IV
%   - Group: DYS only
%
% Analysis 4: Diurnal theta upscaling vs UDysRS Part IV
%   - Group: DYS only
%
% Regions:
%   - SigAll, Motor, Prefrontal
%   - Figure 5 correlations: Prefrontal (F3, F4, Fz)
%
% Statistical analysis:
%   - Covariates: Age and Sex
%   - Pearson correlations on covariate-adjusted residuals
%   - Fisher's z-transformation for 95% confidence intervals
%   - Benjamini-Hochberg FDR correction across 6 LEDD tests
%   - UDysRS correlations reported as exploratory
%
% INPUT:
%   - EEG PSD data: psdIndividual.mat
%   - Clinical metadata: clinical_data.xlsx
%
% OUTPUT:
%   - ThetaMorning_vs_LEDD_results.xlsx
%   - Fig5_All_Correlations.xlsx
%   - Fig5_Correlations_2x4.png
%   - Additional regression plots (.png)
%
% Author: [Giovanni Lombardi]
% Institution: [Institute of Clinical Neuroscience of Southern Switzerland]
% Date: December[2025]
%% ========================================================================

clearvars -except psdIndividual freqs
close all; clc

load('...\psdIndividual.mat');
fileClin   = '...\clinical_data.xlsx';
groupNames = {'EPD','ADV','DYS'};
thetaRange = [4 8];

%% Regions
regions = struct( ...
  'SigAll',        {{'Fz','F4','FC2','C3','C4','CP5','CP1','P7','Pz','O2', ...
                     'F5','F2','F6','FC4','C5','CP3','P5','P2','P6','PO5', ...
                     'PO3','PO4','PO6','TP7','PO8'}}, ...
  'Motor',         {{'C3','C4'}}, ...
  'Prefrontal',    {{'F3','F4','Fz'}} ...
);
names = fieldnames(regions);

%% Load PSD & metadata
[psdArray, labels] = build_PSD_data_from_struct(psdIndividual, [{'CTL'}, groupNames]);
ids      = replace(lower(string(labels.Id)),'ta','x');
sess     = string(labels.Session);
grpAll   = string(labels.Group);
thetaIdx = freqs >= thetaRange(1) & freqs <= thetaRange(2);

chanlocsAll = readlocs('...\ANT_64.ced');
channelsToExclude = [13, 19, 32];
chanlocs = chanlocsAll;
chanlocs(channelsToExclude) = [];
chanLabels = string({chanlocs.labels});

uid = unique(ids);
nUID = numel(uid);

%% Extract log(Theta) and Delta log(Theta)
thetaMorning  = struct(); 
deltaLogTheta = struct();
for r = 1:numel(names)
  thetaMorning.(names{r})  = nan(nUID,1);
  deltaLogTheta.(names{r}) = nan(nUID,1);
end

for i = 1:nUID
  s = uid(i);
  im = find(ids==s & sess=='morning',1);
  ie = find(ids==s & sess=='evening',1);

  if ~isempty(im), Pm = squeeze(psdArray(:,:,im)); end
  if ~isempty(ie), Pe = squeeze(psdArray(:,:,ie)); end

  for r = 1:numel(names)
    nm = names{r};
    idxCh = ismember(chanLabels, regions.(nm));
    if ~isempty(im)
      vm = mean(Pm(idxCh,thetaIdx),'all');
      thetaMorning.(nm)(i) = log(vm);
    end
    if ~isempty(im) && ~isempty(ie)
      ve = mean(Pe(idxCh,thetaIdx),'all');
      deltaLogTheta.(nm)(i) = log(ve) - log(vm);
    end
  end
end

%% Load clinical data
Tclin = readtable(fileClin, 'PreserveVariableNames', true);
Tclin.ID = replace(lower(string(Tclin.ID)), 'ta', 'x');

% Fix AIMS column name if needed
v = Tclin.Properties.VariableNames;
idxA = find(strcmpi(v,'AIMS'), 1);
if isempty(idxA)
    error('AIMS column not found in clinical_data.xlsx');
end
Tclin.Properties.VariableNames{idxA} = 'AIMS';

varsToConvert = {'AIMS','UDysRS_Tot','UPDRS_Tot','UDysRS_parteIV','UPDRS_parteIV','LEDD'};
for vn = varsToConvert
    if iscell(Tclin.(vn{1}))
        Tclin.(vn{1}) = cellfun(@str2double, Tclin.(vn{1}));
    end
end

CT = table(uid, 'VariableNames', {'subject'});
CT.subject        = lower(string(CT.subject));
CT.group          = strings(nUID,1);
CT.age            = nan(nUID,1);
CT.sex            = strings(nUID,1);
CT.LEDD           = nan(nUID,1);
CT.AIMS           = nan(nUID,1);
CT.UPDRS_Tot      = nan(nUID,1);
CT.UPDRS_III      = nan(nUID,1);
CT.UPDRS_parteIV  = nan(nUID,1);
CT.UDysRS_Tot     = nan(nUID,1);
CT.UDysRS_parteIV = nan(nUID,1);

for i = 1:nUID
    s = CT.subject(i);
    im = find(ids==s & sess=='morning',1);
    if ~isempty(im), CT.group(i) = grpAll(im); end
    ic = find(Tclin.ID==s,1);
    if ~isempty(ic)
        CT.age(i)            = Tclin.Age(ic);
        CT.sex(i)            = Tclin.Sex(ic);
        CT.LEDD(i)           = Tclin.LEDD(ic);
        CT.AIMS(i)           = Tclin.AIMS(ic);
        CT.UPDRS_Tot(i)      = Tclin.UPDRS_Tot(ic);
        CT.UPDRS_III(i)      = Tclin.UPDRS_III(ic);
        CT.UPDRS_parteIV(i)  = Tclin.UPDRS_parteIV(ic);
        CT.UDysRS_Tot(i)     = Tclin.UDysRS_Tot(ic);
        CT.UDysRS_parteIV(i) = Tclin.UDysRS_parteIV(ic);
    end
end

% Cleanup
if iscell(CT.LEDD), CT.LEDD = cellfun(@str2double, CT.LEDD); end

keep = ismember(CT.group, groupNames) ...
     & ~isnan(CT.age) ...
     & CT.sex ~= "";
CT = CT(keep, :);

%% 5) Theta_morning vs LEDD (per group & region)
groups   = {'EPD','ADV','DYS'};
regionsL = {'SigAll','Prefrontal'};
alpha    = 0.05;
Z        = norminv(1-alpha/2);

resTab = table('Size',[0 6], ...
    'VariableTypes',{'string','string','double','double','double','double'}, ...
    'VariableNames',{'Group','Region','Rho','P','CI_lower','CI_upper'});

fig = figure('Name','Theta\_morning vs LEDD','Color','w');
tiledlayout(numel(regionsL),numel(groups),'Padding','compact','TileSpacing','compact');

for gi = 1:numel(groups)
    grp = groups{gi};
    sel = CT.group==grp & ~isnan(CT.LEDD);
    n   = sum(sel);
    if n<4, continue; end

    Age  = CT.age(sel);
    SexN = double(categorical(CT.sex(sel))=='M');
    LED  = CT.LEDD(sel);

    for ri = 1:numel(regionsL)
        region = regionsL{ri};
        Yraw   = thetaMorning.(region)(sel);

        valid = ~isnan(Yraw) & ~isnan(LED) & ~isnan(Age) & ~isnan(SexN);
        if sum(valid) < 4, continue; end

        Yraw = Yraw(valid);
        LEDv = LED(valid);
        Agev = Age(valid);
        Sexv = SexN(valid);
        n    = numel(Yraw);

        Tfy = table(Yraw, Agev, Sexv, 'VariableNames',{'Y','Age','Sex'});
        Tfx = table(LEDv, Agev, Sexv, 'VariableNames',{'X','Age','Sex'});

        resY = fitlm(Tfy,'Y~Age+Sex').Residuals.Raw;
        resX = fitlm(Tfx,'X~Age+Sex').Residuals.Raw;

        [rho,p] = corr(resX, resY, 'Type','Pearson');
        z       = atanh(rho);
        se      = 1/sqrt(n-3);
        zCI     = z + Z*[-se, +se];
        rhoCI   = tanh(zCI);

        resTab(end+1,:) = {grp, region, rho, p, rhoCI(1), rhoCI(2)};

        ax = nexttile((ri-1)*numel(groups) + gi);
        scatter(ax, resX, resY, 20, 'filled'); 
        hold(ax,'on');

        tblLM = table(resX, resY, 'VariableNames',{'X','Y'});
        lm    = fitlm(tblLM, 'Y~X');
        xx    = linspace(min(resX),max(resX),100)';
        [yPred,yCI] = predict(lm, table(xx,'VariableNames',{'X'}), 'Alpha', alpha);

        Xb = [xx; flipud(xx)];
        Yb = [yCI(:,1); flipud(yCI(:,2))];
        patch(ax, Xb, Yb, 'r', 'FaceAlpha', .2, 'EdgeColor','none');
        plot(ax, xx, yPred, 'r-', 'LineWidth',1.2);
        title(ax, sprintf('%s – %s\nr=%.2f, p=%.3f', grp, region, rho, p), 'Interpreter', 'none');
    end
end

disp(resTab);
writetable(resTab, 'ThetaMorning_vs_LEDD_results.xlsx');
exportgraphics(gcf, 'ThetaMorning_vs_LEDD_commented.png', 'Resolution', 600);
print(gcf, 'ThetaMorning_vs_LEDD_clean.png', '-dpng', '-r600');

%% 6c) UDysRS_parteIV vs Theta_morning (only DYS)
sel = CT.group=='DYS';

Yall    = CT.UDysRS_parteIV(sel);
thetaS  = thetaMorning.SigAll(sel);
thetaP  = thetaMorning.Prefrontal(sel);
regionsL = {'SigAll','Prefrontal'};
thetav   = {thetaS, thetaP};

Age   = CT.age(sel);
SexN  = double(categorical(CT.sex(sel))=='M');
alpha = 0.05;
Z     = norminv(1-alpha/2);

resTab = table('Size',[0 6], ...
    'VariableTypes',{'string','string','double','double','double','double'}, ...
    'VariableNames',{'Group','Region','Rho','P','CI_lower','CI_upper'});

fig = figure('Name','UDysRS_parteIV vs Theta\_morning (DYS)','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

for ri = 1:2
    region = regionsL{ri};
    Xraw   = thetav{ri};
    Yraw   = Yall;

    valid = ~isnan(Xraw) & ~isnan(Yraw) & ~isnan(Age) & ~isnan(SexN);
    X = Xraw(valid);
    Y = Yraw(valid);
    a = Age(valid);
    s = SexN(valid);

    n = numel(Y);
    if n<3
        warning('Only %d valid for %s', n, region);
        continue;
    end

    T = table(Y, X, a, s, ...
        'VariableNames',{'Y','X','Age','Sex'});

    resY = fitlm(removevars(T,'X'), 'Y ~ Age + Sex').Residuals.Raw;
    resX = fitlm(removevars(T,'Y'), 'X ~ Age + Sex').Residuals.Raw;

    [rho,p] = corr(resX, resY, 'Type','Pearson');

    resTab(end+1,:) = {"DYS", region, rho, p, ...
        tanh(atanh(rho)-Z/sqrt(n-3)), tanh(atanh(rho)+Z/sqrt(n-3))};

    ax = nexttile(ri);
    scatter(ax, resX, resY, 25, 'filled'); 
    hold(ax,'on');

    LM = fitlm(table(resX,resY,'VariableNames',{'X','Y'}),'Y~X');
    xx = linspace(min(resX),max(resX),100)';
    [yPred,yCI] = predict(LM, table(xx,'VariableNames',{'X'}), 'Alpha', alpha);

    patch(ax, [xx; flipud(xx)], [yCI(:,1); flipud(yCI(:,2))], ...
        'b','FaceAlpha',.2,'EdgeColor','none');
    plot(ax, xx, yPred, 'b-','LineWidth',1.5);
    title(ax, sprintf('DYS – %s\nρ=%.2f, p=%.3f', region, rho, p), 'Interpreter', 'none');
end

disp('--- UDysRS_parteIV vs Theta_morning (DYS) ---');
disp(resTab);


%% ========================================================================
% 7) FIGURE 5 - MORNING AND DIURNAL THETA UPSCALING
% ADDITIONAL BLOCK - after revision
% Age- and sex-adjusted residual correlations
% LEDD: 6 tests with BH-FDR correction
% UDysRS Part IV: 2 exploratory correlations
% ========================================================================

% Align theta data with subjects retained in CT
[foundID, idxUID] = ismember(string(CT.subject), string(uid));
assert(all(foundID), 'Subject ID mismatch');

% Frontal region: F3, F4, Fz
thetaM = thetaMorning.Prefrontal(idxUID);
thetaD = deltaLogTheta.Prefrontal(idxUID);

Ydata = {thetaM, thetaD};
outcomes = {'Morning', 'Evening-Morning'};
groups5 = {'EPD','ADV','DYS'};

% Results table
R5 = table('Size',[0 9], ...
    'VariableTypes',{'string','string','string', ...
                     'double','double','double', ...
                     'double','double','double'}, ...
    'VariableNames',{'Clinical','Group','Outcome', ...
                     'N','Rho','P','P_FDR', ...
                     'CI_lower','CI_upper'});

figure('Name','Figure 5 - Theta correlations','Color','w', ...
       'Position',[100 100 1500 700]);

tl = tiledlayout(2,4,'Padding','compact', ...
               'TileSpacing','compact');
axs = gobjects(2,4);

for oi = 1:2

    Yall = Ydata{oi};

    for ci = 1:4

        if ci <= 3
            grp = string(groups5{ci});
            clinical = "LEDD";
            Xall = CT.LEDD;
            colorPlot = [0.9 0.1 0.1];
        else
            grp = "DYS";
            clinical = "UDysRS";
            Xall = CT.UDysRS_parteIV;
            colorPlot = [0.1 0.2 0.9];
        end

        valid = CT.group == grp & ...
                isfinite(Xall) & isfinite(Yall) & ...
                isfinite(CT.age) & CT.sex ~= "";

        n = sum(valid);

        ax = nexttile(tl,(oi-1)*4+ci);
        axs(oi,ci) = ax;

        if n < 5
            warning('Insufficient N: %s %s',grp,outcomes{oi});

            R5(end+1,:) = {clinical,grp,string(outcomes{oi}), ...
                           n,NaN,NaN,NaN,NaN,NaN};
            continue;
        end

        X = double(Xall(valid));
        Y = double(Yall(valid));
        Age = double(CT.age(valid));
        Sex = double(strcmpi(string(CT.sex(valid)),'M'));

        % Same residualization as original script
        Tx = table(X,Age,Sex);
        Ty = table(Y,Age,Sex);

        resX = fitlm(Tx,'X ~ Age + Sex').Residuals.Raw;
        resY = fitlm(Ty,'Y ~ Age + Sex').Residuals.Raw;

        % Pearson correlation of residuals
        [rho,p] = corr(resX,resY,'Type','Pearson');

        % Fisher confidence interval
        z = atanh(rho);
        se = 1/sqrt(n-3);
        zCI = z + norminv(0.975)*[-se,+se];
        rhoCI = tanh(zCI);

        R5(end+1,:) = {clinical,grp,string(outcomes{oi}), ...
                       n,rho,p,NaN,rhoCI(1),rhoCI(2)};

        % Scatterplot and regression
        scatter(ax,resX,resY,25,'filled', ...
                'MarkerFaceColor',[0.1 0.4 0.8]);
        hold(ax,'on');

        LM = fitlm(table(resX,resY, ...
                   'VariableNames',{'X','Y'}),'Y ~ X');

        xx = linspace(min(resX),max(resX),100)';
        [yPred,yCI] = predict(LM, ...
            table(xx,'VariableNames',{'X'}),'Alpha',0.05);

        patch(ax,[xx;flipud(xx)], ...
              [yCI(:,1);flipud(yCI(:,2))], ...
              colorPlot,'FaceAlpha',0.2,'EdgeColor','none');

        plot(ax,xx,yPred,'Color',colorPlot,'LineWidth',1.5);

        if ci <= 3
            xlabel(ax,'Residual LEDD');
        else
            xlabel(ax,'Residual UDysRS Part IV');
        end

        if oi == 1
            ylabel(ax,'Residual Theta MORNING');
        else
            ylabel(ax,'Residual Theta EVENING - MORNING');
        end

    end
end

%% BH-FDR across the 6 LEDD correlations only
idx = find(R5.Clinical=="LEDD" & isfinite(R5.P));
pvals = R5.P(idx);
m = numel(pvals);

if m ~= 6
    warning('Expected 6 LEDD tests; found %d',m);
end

if m > 0
    [pSorted,order] = sort(pvals);

    qSorted = pSorted .* (m./(1:m)');
    qSorted = flipud(cummin(flipud(qSorted)));
    qSorted = min(qSorted,1);

    q = nan(m,1);
    q(order) = qSorted;

    R5.P_FDR(idx) = q;
end

%% Update plot titles with p-values and FDR
for oi = 1:2
    for ci = 1:4

        k = (oi-1)*4+ci;

        if R5.Clinical(k)=="LEDD"
            title(axs(oi,ci),sprintf( ...
                '%s\nr=%.2f, p=%.3f, pFDR=%.3f', ...
                char(R5.Group(k)),R5.Rho(k), ...
                R5.P(k),R5.P_FDR(k)));
        else
            title(axs(oi,ci),sprintf( ...
                'DYS\nr=%.2f, p=%.3f', ...
                R5.Rho(k),R5.P(k)));
        end
    end
end

%% Export
disp('--- Figure 5: all correlations ---');
disp(R5);

writetable(R5,'Fig5_All_Correlations.xlsx');

exportgraphics(gcf,'Fig5_Correlations_2x4.png', ...
               'Resolution',600);
