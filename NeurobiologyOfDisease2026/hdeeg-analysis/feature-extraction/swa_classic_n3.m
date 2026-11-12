%% ========================================================================
%  SWA CLASSIC in N3 per 10 segmenti definiti su NREM (N1+N2+N3)
%
%  IDEA:
%   1) Segmentazione fatta sulla sequenza delle epoche NREM valide
%   2) Per ogni segmento prendo solo le epoche N3
%   3) Concateno i campioni N3
%   4) Calcolo la power spettrale con pwelch
%   5) Estraggo la potenza nella banda SWA
%
%  METRICA:
%   - SWAclassic_N3 = potenza media nella banda 0.5–4 Hz
%
%  NOTE:
%   
%  - Early e Late definiti sui segmenti NREM
%
%  REQUIREMENTS:
%   - EEGLAB
% ========================================================================

clear; clc;

%% ===================== USER SETTINGS ====================================
in_root   = 'G:\EOC\Discinesia\LID-PD-SLEEP\Data\Annals_processed';
score_dir = fullfile(in_root, 'infoscoringpersubject');
out_root  = fullfile(in_root, 'for_reviewers', 'SWAclassic_pwelch');
groups    = {'CTL','DNV','ADV','DYS'};

n_segments         = 10;
early_segments_idx = [2 3 4 5]; % 2345
late_segments_idx  = [6 7 8 9]; % 6789

epoch_len_s = 30;
f_band      = [0.5 4];

exclude_arousals = true;
do_avgref        = false;
do_group_plots   = true;

if ~exist(out_root,'dir'), mkdir(out_root); end

%% ===================== SUMMARY STRUCT ===================================
groupSum = struct();
for g = 1:numel(groups)
    grp = groups{g};
    groupSum.(grp) = struct('subj',{{}}, 'classic_E',[], 'classic_L',[]);
end

%% ===================== LOOP GRUPPI ======================================
for g = 1:numel(groups)

    grp    = groups{g};
    in_dir = fullfile(in_root, grp);
    out_dir = fullfile(out_root, grp);
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    D = dir(fullfile(in_dir,'*.set'));
    if isempty(D), continue; end

    fprintf('\n=== %s: %d file ===\n', grp, numel(D));

    for i = 1:numel(D)

        fname = D(i).name;
        subj  = erase(fname,'.set');

        score_file = fullfile(score_dir, sprintf('%s_%s_scoring.mat', grp, subj));
        out_file   = fullfile(out_dir, sprintf('%s_swaClassic10seg.mat', subj));

        try
            %% -------- LOAD EEG --------
            EEG = pop_loadset('filename', fname, 'filepath', in_dir);
            fs  = EEG.srate;

            if do_avgref
                EEG = pop_reref(EEG, []);
            end

            %% -------- LOAD SCORING --------
            S = load(score_file, 'epochStages', 'epochArousals');
            epochStages   = S.epochStages(:)';
            epochArousals = S.epochArousals(:)';

            %% -------- MATCH EEG / SCORING --------
            if EEG.trials > 1
                nEpochData  = EEG.trials;
                samps_epoch = EEG.pnts;
                data = double(EEG.data);
            else
                samps_epoch = round(fs * epoch_len_s);
                nEpochData  = floor(size(EEG.data,2) / samps_epoch);
                data = double(EEG.data(:,1:nEpochData*samps_epoch));
            end

            nEpochUse = min(nEpochData, numel(epochStages));
            epochStages   = epochStages(1:nEpochUse);
            epochArousals = epochArousals(1:nEpochUse);

            if EEG.trials > 1
                data = data(:,:,1:nEpochUse);
            else
                data = reshape(data, EEG.nbchan, samps_epoch, nEpochUse);
            end

            %% -------- VALID EPOCHS --------
            isN1 = epochStages == 1;
            isN2 = epochStages == 2;
            isN3 = epochStages == 3;

            if exclude_arousals
                isOK = epochArousals == 0;
            else
                isOK = true(size(epochStages));
            end

            valid_NREM = find((isN1 | isN2 | isN3) & isOK);
            valid_N3   = find(isN3 & isOK);

            segEpochs_NREM = segment10_split(valid_NREM, n_segments);

            %% -------- OUTPUT --------
            RES = struct();
            RES.subj = subj;
            RES.group = grp;
            RES.fs = fs;
            RES.f_band = f_band;
            RES.n_segments = n_segments;
            RES.early_segments_idx = early_segments_idx;
            RES.late_segments_idx  = late_segments_idx;
            RES.seg = repmat(struct( ...
                'seg_id',[], ...
                'epochs_NREM',[], ...
                'epochs_N3',[], ...
                'nEpochNREM',0, ...
                'nEpochN3',0, ...
                'SWAclassic_N3', nan(EEG.nbchan,1)), n_segments, 1);

            %% -------- MAIN LOOP --------
            for s = 1:n_segments

                epochsNREM = segEpochs_NREM{s};
                epochsN3   = intersect(epochsNREM, valid_N3, 'stable');

                RES.seg(s).seg_id      = s;
                RES.seg(s).epochs_NREM = epochsNREM;
                RES.seg(s).epochs_N3   = epochsN3;
                RES.seg(s).nEpochNREM  = numel(epochsNREM);
                RES.seg(s).nEpochN3    = numel(epochsN3);

                if isempty(epochsN3), continue; end

                for ch = 1:EEG.nbchan
                    Xn3  = squeeze(data(ch,:,epochsN3));   % [samples x epochs]
                    xvec = Xn3(:);

                    [pxx,f] = pwelch(xvec, [], [], [], fs);
                    RES.seg(s).SWAclassic_N3(ch) = bandpower_from_psd(pxx, f, f_band);
                end
            end

            %% -------- SAVE --------
            save(out_file, 'RES', '-v7');

            %% -------- GROUP SUMMARY --------
            mClassic = seg_means_over_channels(RES);

            groupSum.(grp).subj{end+1,1}   = subj;
            groupSum.(grp).classic_E(end+1,1) = mean(mClassic(early_segments_idx), 'omitnan');
            groupSum.(grp).classic_L(end+1,1) = mean(mClassic(late_segments_idx),  'omitnan');

            fprintf(' Saved: %s\n', out_file);

        catch ME
            warning(' %s %s: %s', grp, subj, ME.message);
        end
    end

    %% -------- GROUP PLOT --------
    if do_group_plots && ~isempty(groupSum.(grp).classic_E)
        GS = groupSum.(grp);
        figure('Color','w','Name',sprintf('%s | SWA classic Early vs Late', grp));
        paired_plot(GS.classic_E, GS.classic_L, 'SWA classic (0.5–4 Hz)');
        sgtitle(sprintf('%s | Early vs Late', grp));
    end
end

%% ===================== LOCAL FUNCTIONS ==================================
function seg = segment10_split(valid_epochs, nseg)
    seg = cell(1,nseg);
    if isempty(valid_epochs)
        for s = 1:nseg, seg{s} = []; end
        return;
    end
    K  = numel(valid_epochs);
    bx = round(linspace(1, K+1, nseg+1));
    for s = 1:nseg
        if bx(s) <= bx(s+1)-1
            seg{s} = valid_epochs(bx(s):bx(s+1)-1);
        else
            seg{s} = [];
        end
    end
end

function bp = bandpower_from_psd(pxx, f, band)
    idx = f >= band(1) & f <= band(2);
    bp  = trapz(f(idx), pxx(idx));
end

function mClassic = seg_means_over_channels(RES)
    nseg = numel(RES.seg);
    mClassic = nan(nseg,1);
    for s = 1:nseg
        mClassic(s) = mean(RES.seg(s).SWAclassic_N3, 'omitnan');
    end
end

function paired_plot(E, L, ylab)
    hold on;
    for i = 1:numel(E)
        plot([1 2], [E(i) L(i)], '-o');
    end
    xlim([0.7 2.3]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Early','Late'});
    ylabel(ylab);
    grid on;
    title(ylab);
end