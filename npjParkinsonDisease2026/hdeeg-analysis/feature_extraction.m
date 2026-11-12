%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEG Power Spectral Density (PSD) Computation Script
%
% Description:
% This script processes preprocessed EEG `.set` files organized by group 
% (CTL, EPD, ADV, DYS), calculates power spectral density (PSD) for each 
% subject and session (morning/evening), applies a notch interpolation 
% around 50 Hz, and stores results in a structured variable for further 
% analysis. Channel 32 (EOG) is removed before computation.
%
% Output:
% - PSD values per subject/channel stored in `psdIndividual.mat`
%
% Requirements:
% - EEGLAB installed and initialized
% - Input `.set` files stored in per-group folders
% Author: [Giovanni Lombardi]
% Institution: [Institute of Clinical Neuroscience of Southern Switzerland]
% Date: December[2025]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial configuration
clear; clc;
eeglabPath = '...eeglab_current\eeglab2024.2';
addpath(eeglabPath);
eeglab;

% Define paths
% dataPath = '...data_preprocessed_BSL';
dataPath = '...._full_eeg_downsamp_data_fs_250_meanRef';
outputPath = fullfile(dataPath, 'dati_PSD');  % Output folder

% Create output folder if it doesn't exist
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

% Define subject groups
groups = {'CTL', 'EPD', 'ADV', 'DYS'};

% Analysis parameters
window = 512;                        % Window length (samples)
noverlap = round(window / 2);       % 50% overlap
nfft = window;                      % FFT points
notchBand = [49.5, 50.5];           % Notch band around 50 Hz

% Initialize structure to store results
psdIndividual = struct();

%% Process each group
for g = 1:length(groups)
    group = groups{g};
    disp(['Processing group: ', group]);
    
    groupOutputFolder = fullfile(outputPath, group);
    if ~exist(groupOutputFolder, 'dir')
        mkdir(groupOutputFolder);
    end

    groupFolder = fullfile(dataPath, group);
    allFiles = dir(fullfile(groupFolder, '*.set'));
    if isempty(allFiles)
        warning(['No files found for group ', group]);
        continue;
    end
    
    for f = 1:length(allFiles)
        filePath = fullfile(groupFolder, allFiles(f).name);
        
        % Determine session (morning or evening)
        session = 'morning';
        if contains(allFiles(f).name, '_s1_')
            session = 'evening';
        end
        
        % Load EEG file
        disp(['Loading file: ', filePath]);
        EEG = pop_loadset(filePath);

        % Remove channel 32 (EOG)
        disp(['Removing channel 32 (EOG) - ', session, '...']);
        EEG = pop_select(EEG, 'nochannel', {'32'});
        EEG = eeg_checkset(EEG);
        
        % Compute PSD per channel
        psdTemp = [];
        for ch = 1:EEG.nbchan
            [pxx, freqs] = pwelch(EEG.data(ch, :), window, noverlap, nfft, EEG.srate);
            
            % Identify notch region
            notchIdx = (freqs >= notchBand(1) & freqs <= notchBand(2));
            freqsBefore = freqs(find(notchIdx, 1) - 1);
            freqsAfter = freqs(find(notchIdx, 1, 'last') + 1);
        
            % Interpolate notch region
            pxxInterp = pxx;
            pxxInterp(notchIdx) = interp1([freqsBefore, freqsAfter], ...
                                          [pxx(find(notchIdx, 1) - 1), pxx(find(notchIdx, 1, 'last') + 1)], ...
                                          freqs(notchIdx), 'nearest');
        
            % Append interpolated PSD
            psdTemp = [psdTemp; pxxInterp'];
        end
        
        % Save results into structured variable
        subjectID = matlab.lang.makeValidName(allFiles(f).name(1:end-4)); % Clean file name
        disp(subjectID);

        if ~isfield(psdIndividual, group)
            psdIndividual.(group) = struct();
        end
        if ~isfield(psdIndividual.(group), session)
            psdIndividual.(group).(session) = struct();
        end
        psdIndividual.(group).(session).(subjectID) = psdTemp;
    end
end

% Final save
disp('Saving individual PSD results...');
save(fullfile(outputPath, 'psdIndividual.mat'), 'psdIndividual', 'freqs', '-v7.3');
disp('PSD computation completed.');
