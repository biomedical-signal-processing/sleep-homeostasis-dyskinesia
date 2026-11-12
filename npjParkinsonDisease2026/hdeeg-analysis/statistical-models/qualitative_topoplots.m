%% ========================================================================
%  EEG ANALYSIS – Theta Power Topoplots
%  - Morning, Evening and Evening-Morning qualitative comparison across 4 groups (CTL, EPD, ADV, DYS)
%  - Theta band analyzed: 4–8 Hz
%  - Output: log10 power topoplots and difference maps (evening - morning)
%  - Results are saved in the "report" subfolder
% Author: [Giovanni Lombardi]
% Institution: [Institute of Clinical Neuroscience of Southern Switzerland]
% Date: December[2025]
%  ========================================================================

clear; clc;

% === Define paths and parameters ===
dataPath = '...brain_full_fs_250_meanRef';
groups   = {'CTL', 'EPD', 'ADV', 'DYS'};

% Define theta ands to analyze
% eg:[4 8], [4 6], [5 7], [6 8]
thetaBand = [4 8];
channelsToRemove = [13, 19, 32];  % unused channels

% === Load PSD data for all groups ===
disp('Loading saved PSD data...');
psdData = struct();
chanlocs = [];

for g = 1:length(groups)
    group = groups{g};
    groupPath = fullfile(dataPath, group);

    load(fullfile(groupPath, 'psdGroup_05_100Hz_w512_spline.mat'), ...
        'psdGroup', 'freqs', 'EEG_morning');

    psdData.(group).morning = psdGroup.(group).morning.mean;
    psdData.(group).evening = psdGroup.(group).evening.mean;

    psdData.(group).morning(channelsToRemove, :) = [];
    psdData.(group).evening(channelsToRemove, :) = [];

    if isempty(chanlocs)
        chanlocs = EEG_morning.chanlocs;
        chanlocs(channelsToRemove) = [];
    end
end

%% === Compute and plot topomaps for theta band ===
fprintf('\nGenerating theta band topoplots: %.1f–%.1f Hz\n', thetaBand(1), thetaBand(2));

% Frequency index for theta band
idxFreq = freqs >= thetaBand(1) & freqs <= thetaBand(2);

% Compute log10 power for morning, evening, and their difference
powerBands = struct();
globalDiffLog = [];

for g = 1:length(groups)
    group = groups{g};

    morning = mean(psdData.(group).morning(:, idxFreq), 2);
    evening = mean(psdData.(group).evening(:, idxFreq), 2);

    powerBands.(group).morning = log10(morning + eps);
    powerBands.(group).evening = log10(evening + eps);
    powerBands.(group).diff    = powerBands.(group).evening - powerBands.(group).morning;

    globalDiffLog = [globalDiffLog; powerBands.(group).diff];
end

% Define common color scale for difference maps
minScale = min(globalDiffLog);
maxScale = max(globalDiffLog);

% === Create figure with topoplots ===
nRows = length(groups);
nCols = 3;  % morning | evening | difference
figure('Name', sprintf('Theta Band %.1f–%.1f Hz', thetaBand(1), thetaBand(2)), ...
       'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.85]);

t = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

for g = 1:nRows
    group = groups{g};

    morningData = powerBands.(group).morning;
    eveningData = powerBands.(group).evening;
    clim_common = [min([morningData; eveningData]), max([morningData; eveningData])];

    for c = 1:nCols
        ax = nexttile((g-1)*nCols + c);

        switch c
            case 1
                data = morningData;
                clim = clim_common;
            case 2
                data = eveningData;
                clim = clim_common;
            case 3
                data = powerBands.(group).diff;
                clim = [minScale, maxScale];
        end

        topoplot(data, chanlocs, 'maplimits', clim, ...
            'electrodes', 'on', 'style', 'map', 'plotrad', 0.5);
        colorbar;
        axis off;
    end
end

%% === Save figure ===
outputFolder = fullfile(dataPath, 'report');
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

fname = sprintf('theta_%d-%dHz_clean', round(thetaBand(1)), round(thetaBand(2)));
saveas(gcf, fullfile(outputFolder, [fname, '.png']));
savefig(gcf, fullfile(outputFolder, [fname, '.fig']));
print(gcf, fullfile(outputFolder, [fname, '.tiff']), '-dtiffn', '-r600');

fprintf(' Figure saved: %s\n', fullfile(outputFolder, fname));
