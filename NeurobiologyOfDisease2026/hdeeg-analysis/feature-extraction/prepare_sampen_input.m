%% ====== SETTINGS ======
base_dir = '...\Annals_processed\ADV';                    % folder containing the .set files
out_dir  = '...\entropy_input_data\ADV';                 % output folder for minimal .mat files

% Optional: create output folder if it does not exist
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% ====== EEGLAB on path (if needed) ======
% addpath('...\eeglab2024');
% eeglab nogui;

%% ====== Loop over .set files and save one subject at a time ======
files = dir(fullfile(base_dir, '*.set'));
assert(~isempty(files), 'No .set files found in: %s', base_dir);

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    [~, id_noext] = fileparts(fpath);   % filename without extension

    try
        EEG = pop_loadset('filename', files(k).name, 'filepath', files(k).folder);

        % --- Minimal structure required for downstream analyses ---
        S = struct();
        S.data = EEG.data;    % [channels x samples]
        S.fs   = EEG.srate;   % sampling frequency in Hz

        % Optional metadata
        % S.nchan = size(EEG.data, 1);
        % S.nsamp = size(EEG.data, 2);

        % Output filename for this subject
        out_file = fullfile(out_dir, sprintf('%s_eeg_minimal.mat', id_noext));

        % Save only the minimal structure S for this subject
        if numel(EEG.data) > 2e9
            save(out_file, 'S', '-v7.3');
        else
            save(out_file, 'S');
        end

        fprintf('✔ Saved %s (fs = %.2f Hz)\n', out_file, S.fs);

    catch ME
        warning('Error processing %s: %s', files(k).name, ME.message);
    end

    % Clear memory between iterations
    clear EEG S
end