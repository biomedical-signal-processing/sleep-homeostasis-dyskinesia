function [psdArray, labels_study] = build_PSD_data_from_struct(psdIndividual, groupNames)

groupNames = {'CTL', 'EPD', 'ADV', 'DYS'};
% Channels to exclude
channelsToExclude = [13, 19, 32];

psdList = {};
labelsList = {};

for g = 1:length(groupNames)
    group = groupNames{g};
    
    for session = ["morning", "evening"]
        subjStruct = psdIndividual.(group).(session);
        subjNames = fieldnames(subjStruct);
        
        for s = 1:length(subjNames)
            subjIdFull = subjNames{s};
            data = subjStruct.(subjIdFull); % [channels x freqs]

            % Exclude specified channels
            data(channelsToExclude, :) = [];

            % Extract subject ID from the field name
            tokens = regexp(subjIdFull, 'x(\d+)', 'tokens');
            if ~isempty(tokens)
                subjId = ['x' tokens{1}{1}];
            else
                subjId = subjIdFull;
            end

            % Save PSD data and metadata
            psdList{end+1} = data; % [ch x freq]
            labelsList{end+1,1} = subjId;
            labelsList{end,2} = char(session);
            labelsList{end,3} = group;
        end
    end
end

% Concatenate PSD array: [channels x freqs x subjects]
psdArray = cat(3, psdList{:});

% Create metadata table
labels_study = table(...
    categorical(labelsList(:,1)), ...
    categorical(labelsList(:,2)), ...
    categorical(labelsList(:,3)), ...
    'VariableNames', {'Id', 'Session', 'Group'});

end
