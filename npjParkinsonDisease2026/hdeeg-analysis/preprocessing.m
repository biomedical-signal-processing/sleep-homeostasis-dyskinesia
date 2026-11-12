%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUNLOG THETA ANALYSIS PIPELINE
%
% Description:
% This script performs full preprocessing of raw EEG data for theta-band
% analyses. It is intended for use in clinical EEG datasets organized by
% subject, condition, and group. The preprocessing includes:
%
% 1. Initialization and setup of EEGLAB and dependencies
% 2. File import and structure setup
% 3. Basic filtering (notch, high-pass, low-pass)
% 4. Noisy channel detection and rejection
% 5. Artifact Subspace Reconstruction (ASR) to remove transient noise
% 6. ICA decomposition and component classification (ICLabel)
% 7. Semi-automatic ICA component rejection (visual + ICLabel)
% 8. Visual inspection and removal of noisy segments
% 9. Final saving of the cleaned EEG dataset
%
% Output:
% - Cleaned `.set` files stored per group and condition
% - Intermediate files with ASR and ICA results also saved
%
% Requirements:
% - EEGLAB installed with required plugins:
%   ICLabel, clean_rawdata, viewprops, FieldTrip
% - Proper organization of input folders per subject and group
% - Manual inspection via GUI is required in steps 6, 7, and 8
%
% Author: [Giovanni Lombardi]
% Institution: [Institute of Clinical Neuroscience of Southern Switzerland]
% Date: December[2025]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% --- 1. Initial configuration --- run once to initialize
% Start EEGLAB, then comment this out
% Change to the main folder
unitdisk = 'D';
mainpath = ':\ThetaProject\';
cd([unitdisk mainpath 'EEGDATA\']);

% Add paths to required scripts
restoredefaultpath;
addpath([unitdisk mainpath]);
addpath([unitdisk mainpath 'SCRIPTS\']);
addpath([unitdisk mainpath 'SCRIPTS\ICLabel-master\']);
addpath([unitdisk mainpath 'SCRIPTS\viewprops-master\']);
addpath([unitdisk mainpath 'SCRIPTS\fieldtrip\']);
addpath(genpath([unitdisk mainpath 'SCRIPTS\eeglab2024.2\']));
addpath(genpath([unitdisk mainpath 'SCRIPTS\ANTeepimport1.13\']));
addpath(genpath([unitdisk mainpath 'SCRIPTS\libeep-3.3.177-matlab\']));
addpath(genpath([unitdisk mainpath 'SCRIPTS\csc-eeg-tools-develop\']));
eeglab;

%clear; clc; close all;

%% --- 1. Initial configuration ---
% Change to the main folder
unitdisk = 'D';
mainpath = ':\ThetaProject\';
cd([unitdisk mainpath 'EEGDATA\']);

%% --- 2. Initial parameters ---
ID = 'P065';
subjectID = '065';  % Subject number
condition = '*s1_oc';  % Condition (e.g., morning eyes open)
group = 'EPD';  % Subject group: ctl, epd, adv, dys
dataPath = '....data_processed_BSL';
preprocessedFileName = [subjectID '_' condition '_preprocessedBSL.set'];

% Full path to save group data
groupFolder = fullfile(dataPath, group);

% Create directory if it doesn't exist
if ~exist(groupFolder, 'dir')
    mkdir(groupFolder);
    disp(['Created folder for group ' group ': ' groupFolder]);
else
    disp(['Folder already exists: ' groupFolder]);
end

%% --- 3. File paths and electrode montage ---
eegpath = [unitdisk mainpath 'EEGDATA\' ID '\BSL\'];
%eegpath = [unitdisk mainpath 'EEGDATA\P002\BSL\'];

eegfiles = dir([eegpath 'original\' condition '*.cnt']);
eegname = eegfiles.name;
cd(eegfiles.folder);
plotpath = fullfile(eegpath, 'plots');
mkdir(plotpath);

% Read electrode montage
chanlocs = readlocs([unitdisk mainpath 'SCRIPTS\ANT_64.ced']);

%% --- 4. EEG file import ---
disp('Importing EEG file...');
ep_length = 30;
ft_defaults;
properties = ft_read_header(fullfile(eegfiles.folder, eegname)); 
endrec = floor(properties.nSamples / properties.Fs / ep_length);
channels = ft_read_data(eegname, 'begsample', 1, ...
             'endsample', endrec * ep_length * properties.Fs);
if isempty(channels)
    error('The imported file is not valid.');
end

% Create EEG set file
disp('Creating EEG set file...');
EEG = eeg_emptyset;
EEG.data = channels(1:64, :);
EEG.filename = [ID '_' condition(2:end)];
EEG.filepath = eegpath;
EEG.setname = eegname(1:end-4);
EEG.srate = properties.Fs;
EEG.pnts = size(EEG.data, 2);
EEG.nbchan = size(EEG.data, 1);
EEG.chanlocs = chanlocs;
EEG.times = (1:EEG.pnts) / EEG.srate;
EEG.ref = 'monopolar';
EEG.trials = 1;
EEG.xmin = 0;
EEG.xmax = properties.nSamples;
EEG.labels = properties.label;
EEG = eeg_checkset(EEG);

% Save original EEG file
disp('Saving original EEG file...');
EEG = pop_saveset(EEG, 'filename', EEG.filename, 'filepath', EEG.filepath);


%% --- 5. Signal preprocessing ---
disp('Signal preprocessing...');

% Apply filtering
disp('Applying filters...');
EEG = pop_eegfiltnew(EEG, 'locutoff', 49.5, 'hicutoff', 50.5, 'revfilt', 1);  % Notch filter at 50 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);  % High-pass filter at 0.5 Hz
%EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);  % Optional low-pass filter at 40 Hz

% Check for NaN or Inf after filtering
if any(isnan(EEG.data(:))) || any(isinf(EEG.data(:)))
    error('Filtered signal contains NaN or Inf values.');
end

EEG = eeg_checkset(EEG);

%% --- 6. Detection and handling of noisy channels ---
disp('Detecting noisy channels...');
EEG = csc_eeg_plotter_PD(EEG);  % Manual inspection of noisy channels
EEG.badchannels_manual = [32 EEG.hidden_channels];  % Mark bad channels (e.g., reference + hidden)
EEG.goodchannels_manual = setdiff(1:EEG.nbchan, EEG.badchannels_manual);  % Get valid channels

% Sanity check on channel indices
if max(EEG.goodchannels_manual) > EEG.nbchan
    error('Indices in EEG.goodchannels_manual exceed the number of available channels.');
end
disp(['Updated number of good channels: ', num2str(length(EEG.goodchannels_manual))]);

disp('Press any key to apply 100 Hz low-pass filter...');
pause

% Select only good channels
EEG_good = pop_select(EEG, 'channel', EEG.goodchannels_manual);

% Apply low-pass filter at 100 Hz
disp('Applying 100 Hz low-pass filter...');
EEG_good = pop_eegfiltnew(EEG_good, 'hicutoff', 100);

% Check for NaN or Inf after filtering
if any(isnan(EEG_good.data(:))) || any(isinf(EEG_good.data(:)))
    error('Filtered signal contains NaN or Inf values.');
end

% Final consistency check
EEG_good = eeg_checkset(EEG_good);

%% Artifact removal with ASR
% Parameter configuration for clean_artifacts
% flatlineCriterion = [];   % No flatline threshold
% channelCriterion  = 0.8;  % Remove channels with variance > 80% compared to others
% burstCriterion    = 5;    % Remove bursts with a threshold of 5
% windowCriterion   = [];   % Use default values
% maxBadChannels    = [];   % Use default values

disp('Artifact removal using ASR...');
disp('Removing transient noise bursts with ASR on good channels...');

%EEG_good = pop_select(EEG, 'channel', EEG.goodchannels_manual);
EEG_good = clean_rawdata(EEG_good, -1, -1, -1, -1, 5, -1);

% Check for NaN or Inf values after ASR
if any(isnan(EEG_good.data(:))) || any(isinf(EEG_good.data(:)))
    error('EEG data contains NaN or Inf after ASR.');
end

% Interpolate bad channels back into the dataset
EEG = pop_interp(EEG_good, EEG.chanlocs, 'spherical');

% Save dataset after ASR
EEG = pop_saveset(EEG, 'filename', [EEG.filename '_ASR'], 'filepath', EEG.filepath);

disp('Press any key to proceed with ICA...');
pause

%% --- 7. ICA and component classification ---
disp('Running ICA...');
EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', EEG.goodchannels_manual);

disp('Classifying ICA components using ICLabel...');
EEG = pop_iclabel(EEG, 'default');

disp('Opening interactive GUI for ICA component review...');
% Remove any previously stored hidden channels
if isfield(EEG, 'hidden_channels')
    EEG = rmfield(EEG, 'hidden_channels');
end

pop_viewprops(EEG, 0); % Open GUI for detailed ICA component inspection
EEG = csc_eeg_plotter_PD(EEG);
EEG.components_rej1 = EEG.hidden_channels;
disp('Visually selected components:');
disp(EEG.components_rej1);
disp('Press any key to continue after component review.');
pause;

% Identify ICA components for rejection based on ICLabel
disp('Classifying ICA components...');
classifications = EEG.etc.ic_classification.ICLabel.classifications;
badICs_auto = unique([
    find(classifications(:, 3) > 0.8);  % Eye artifacts
    find(classifications(:, 2) > 0.8);  % Muscle artifacts
    find(classifications(:, 4) > 0.7);  % Heart activity
    find(classifications(:, 6) > 0.8);  % Channel noise
    find(classifications(:, 7) > 0.9)   % Other
]);
badICs_combined = unique([EEG.components_rej1, badICs_auto']);
disp('Total components selected for removal:');
disp(badICs_combined);

EEG = pop_subcomp(EEG, badICs_combined, 0);

% Save cleaned dataset after ICA
EEG = pop_saveset(EEG, 'filename', [EEG.filename '_ICA_cleaned'], 'filepath', EEG.filepath);

disp('Press any key to open the GUI for visual inspection...');
pause

%% --- 8. Visual inspection and removal of noisy segments after ICA ---

% Open GUI for visual inspection of the EEG signal
pop_eegplot(EEG, 1, 1, 1); 

% Remove segments manually marked during GUI inspection
if ~isempty(EEG.reject.rejmanual)
    disp('Removing manually marked segments...');
    EEG = eeg_eegrej(EEG, EEG.reject.rejmanual);  % Remove selected segments
    EEG = eeg_checkset(EEG);
    disp('Noisy segments removed.');
else
    disp('No manually marked segments to remove.');
end

disp('Press any key to proceed with final saving...');
pause

%% --- 9. Final saving ---
disp('Final saving...');
condition_clean = strrep(condition, '*', '_');  % Replace '*' with '_'
EEG = pop_saveset(EEG, 'filename', [subjectID '_' condition_clean '_preprocessedBSL.set'], ...
                  'filepath', groupFolder);

% % Ensure the EEG data is not saved as a separate .fdt file
% EEG = eeg_checkset(EEG, 'makeur'); % Forces the data into memory (non-memory mapped)
% 
% % Save the .set file as a standalone file
% EEG = pop_saveset(EEG, 'filename', [subjectID '_' condition_clean '_preprocessedBSL.set'], ...
                  %'filepath', groupFolder, 'savemode', 'onefile');

disp('Preprocessing completed.');

