%%
% Preprocessing for EEG data:
% - channel reduction from 256 to ~150 channels
% - save subject-level structure (data + fs) for downstream cluster/entropy analyses

clear all; clc

%% === 1. Load .set file ===
[filename, pathname] = uigetfile('*.set');
EEG = pop_loadset('filename', filename, 'filepath', pathname);

%% === Define ROI channel indices (1-based) ===
fp = sort([27,33,34,38,39,47,48,26,20,19,12,11,3,2,222]);
f  = sort([16,22,23,24,28,29,30,35,36,40,41,42,49,50,21,15,7,14,6,207,13,5,215,4,224,223,214,206,213,205]);
c  = sort([9,17,43,44,45,51,52,53,57,58,59,60,64,65,66,71,72,8,81,186,198,197,185,132,196,184,144,204,195,183,155,194,182,164,181,173]);
t  = sort([55,56,62,63,69,70,74,75,84,85,96,221,212,211,203,202,193,192,180,179,171,170]);
p  = sort([76,77,78,79,80,86,87,88,89,97,98,99,100,110,90,101,119,172,163,154,143,131,162,153,142,130,161,152,141,129,128]);
o  = sort([107,108,109,116,117,118,125,126,160,151,140,150,139,127,138]);

roi_map = struct('fp',{fp}, 'f',{f}, 'c',{c}, 't',{t}, 'p',{p}, 'o',{o});
subset_idx = unique([fp, f, c, t, p, o]);

% Safety checks
subset_idx = subset_idx(subset_idx >= 1 & subset_idx <= EEG.nbchan);
assert(~isempty(subset_idx), 'Channel subset is empty with respect to nbchan=%d.', EEG.nbchan);

% Save full original layout for possible future interpolation
origChanlocs_full = EEG.chanlocs;
nbchan_full       = EEG.nbchan;

% Immediately select the channel subset:
% from this point onward, filtering / ASR / ICA run on ~150 channels
EEG = pop_select(EEG, 'channel', subset_idx);
EEG = eeg_checkset(EEG);
fprintf('➡️  Working on channel subset: %d/%d\n', EEG.nbchan, nbchan_full);

%% === 2. Add scoring information (sleep stages + arousals) ===
% Extract subject ID (example: 'PD009')
[~, name_no_ext, ~] = fileparts(filename);
subject_id = regexp(name_no_ext, 'PD\d{3}', 'match', 'once');

if ~isempty(subject_id)

    % === Scoring folders ===
    staging_folder = '...\Staging 128Hz';
    arousal_folder = '...\Final Input for A7';

    % === 2.1 Load stages and arousals ===
    stages_file   = fullfile(staging_folder, ['Stages' subject_id '.mat']);
    arousals_file = fullfile(arousal_folder, ['Arousals' subject_id '.mat']);

    if isfile(stages_file)
        S = load(stages_file);
        A = load(arousals_file);

        if isfield(S, 'stages')
            EEG.stages   = S.stages;
            EEG.arousals = A.arousals;
            fprintf('✅ Stages loaded for %s\n', subject_id);
            fprintf('✅ Arousals loaded for %s\n', subject_id);
        else
            warning('File %s does not contain variable "stages".', stages_file);
        end
    else
        warning('STAGES file not found for %s', subject_id);
    end
else
    warning('Unable to determine subject ID from "%s"', filename);
end

%% === 3. Band-pass filtering (0.5 - 40 Hz) ===
disp('Starting filtering...');
EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);
EEG = pop_eegfiltnew(EEG, 'hicutoff', 40);

% Optional notch filter
EEG = pop_eegfiltnew(EEG, 'locutoff', 49, 'hicutoff', 51, 'revfilt', 1);
EEG = eeg_checkset(EEG);

disp('Filtering completed.');

%% === 4. Re-reference ===
EEG = pop_reref(EEG, []);
disp('Re-referencing completed.');

%% === 5. Identify and remove noisy channels ===
disp('Identifying noisy channels...');
EEG = csc_eeg_plotter_PD(EEG);  % manual channel selection

% Channels marked for removal
EEG.badchannels_manual  = unique([EEG.hidden_channels(:)']);
EEG.goodchannels_manual = setdiff(1:EEG.nbchan, EEG.badchannels_manual);

disp(['Number of good channels: ', num2str(length(EEG.goodchannels_manual))]);

% Save original layout for interpolation
origChanlocs = EEG.chanlocs;

% Keep only good channels
EEG_good = pop_select(EEG, 'channel', EEG.goodchannels_manual);
EEG_good = eeg_checkset(EEG_good);

%% === 6. Epoching into 30-second windows ===
epochLen  = 30;
sampEpoch = epochLen * EEG.srate;
nEpochs   = floor(EEG.pnts / sampEpoch);

%% === 7. ASR for burst artifact removal (currently disabled) ===
% disp('Removing artifacts with ASR...');
% EEG_good = clean_rawdata(EEG_good, -1, -1, -1, -1, 5, -1);
%
% % Quality check
% if any(isnan(EEG_good.data(:))) || any(isinf(EEG_good.data(:)))
%     error(' EEG data contain NaN or Inf values after ASR.');
% end
% disp('ASR completed.');

%% === 8. Interpolate removed channels ===
disp(' Interpolating removed channels...');
EEG = pop_interp(EEG_good, origChanlocs, 'spherical');

%% === 9. Create dummy events for epoching ===
EEG.event = [];
for i = 1:nEpochs
    EEG.event(i).type     = 'X';
    EEG.event(i).latency  = (i-1)*sampEpoch + 1;
    EEG.event(i).duration = sampEpoch;
    EEG.event(i).urevent  = i;
end

EEG = eeg_checkset(EEG, 'eventconsistency');

% Epoching (no baseline removal)
EEG = pop_epoch(EEG, {'X'}, [0 epochLen], 'epochinfo', 'yes');
disp(['Generated ', num2str(EEG.trials), ' epochs of ', num2str(epochLen), ' s.']);

%% === 10. Identify artifact epochs and arousal epochs ===
EEG = csc_eeg_plotter(EEG);   % manual epoch marking
badEp = find(EEG.csc_marked_trials);

% Build scoring aligned with epochs
L = sampEpoch;
S = double(EEG.stages(:));
A = EEG.arousals(:);
N = EEG.trials;

if length(S) < N*L || length(A) < N*L
    warning('stages/arousals shorter than epoched signal. Truncating to %d epochs.', ...
        floor(min([length(S), length(A)]) / L));
    N = floor(min([length(S), length(A)]) / L);
end

S = S(1:N*L);
A = A(1:N*L);

epochStages   = arrayfun(@(k) mode(S((k-1)*L+1 : k*L)), 1:N);
epochArousals = arrayfun(@(k) any(A((k-1)*L+1 : k*L)), 1:N);

arousalEp = find(epochArousals);  % epochs containing arousals

% Merge epochs to reject
badToRemove = union(badEp, arousalEp);
disp(['Epochs to remove (artifact + arousals): ', num2str(length(badToRemove))]);

% Reject epochs from signal
EEG = pop_rejepoch(EEG, badToRemove, 0);

% Remove corresponding scoring entries
epochStages(badToRemove)   = [];
epochArousals(badToRemove) = [];

% Save updated scoring
EEG.swa_scoring.epochStages   = epochStages;
EEG.swa_scoring.epochArousals = epochArousals;
EEG.swa_scoring.badEpochs     = badToRemove;
EEG.swa_scoring.epochCount    = EEG.trials;

disp('Epochs removed and scoring updated.');

%% === 11. Save clean epoched file before ICA ===
EEG.filepath = '...\Annals_processed\ADV';   % change path as needed
EEG = pop_saveset(EEG, ...
    'filename', [EEG.filename 'clean_epoched_beforeICA_'], ...
    'filepath', EEG.filepath);

%% === 12. ICA and component classification ===
disp('Running ICA...');

% ICA on original good channels only
if isfield(EEG, 'goodchannels_manual')
    EEG = pop_runica(EEG, 'icatype', 'runica', 'chanind', EEG.goodchannels_manual);
else
    warning('goodchannels_manual not found. ICA will run on all channels.');
    EEG = pop_runica(EEG, 'icatype', 'runica');
end

%% === 13. ICLabel component classification ===
disp('Classifying ICA components with ICLabel...');
EEG = pop_iclabel(EEG, 'default');

%% === 14. Interactive ICA component inspection ===
disp(' Opening interactive interface for ICA component review...');

% Remove any previous hidden_channels field
if isfield(EEG, 'hidden_channels')
    EEG = rmfield(EEG, 'hidden_channels');
end

pop_viewprops(EEG, 0);      % GUI for component inspection
EEG = csc_eeg_plotter_PD(EEG);  % manual component selection

% Save manually rejected components
EEG.components_rej1 = EEG.hidden_channels;

disp('Components manually marked for rejection:');
disp(EEG.components_rej1);

disp('Press any key to continue after component review.');
pause;

%% === 15. Automatic rejection of artifact-related ICA components ===
disp('Automatic ICA component classification...');

classifications = EEG.etc.ic_classification.ICLabel.classifications;

badICs_auto = unique([
    find(classifications(:, 3) > 0.8);  % Eye
    find(classifications(:, 2) > 0.8);  % Muscle
    find(classifications(:, 4) > 0.8);  % Heart
    find(classifications(:, 6) > 0.9);  % Channel noise
    find(classifications(:, 7) > 0.9)   % Other
]);

% Combine manual + automatic rejection
badICs_combined = unique([EEG.components_rej1, badICs_auto]);

disp('Total ICA components to remove:');
disp(badICs_combined);

% Remove selected components
EEG = pop_subcomp(EEG, badICs_combined, 0);

%% === 16. Visual inspection after ICA ===
disp('Press any key. Opening GUI for post-ICA visual inspection...');
pause;
csc_eeg_plotter(EEG);

%% === 17. Save final cleaned file ===
EEG.filepath = '...\Annals_processed\CTL';   % change path as needed
EEG = pop_saveset(EEG, ...
    'filename', [EEG.filename 'epoched_ICA2_cleaned'], ...
    'filepath', EEG.filepath);