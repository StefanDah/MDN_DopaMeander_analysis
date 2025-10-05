function concat_analysis = concatenate_analysis(analysis)
% GROUPANDCONCATANALYSIS Groups trials and concatenates fields from analysis struct
%
%   concat_analysis = GROUPANDCONCATANALYSIS(analysis)
%
%   Inputs:
%       analysis - Struct array with fields including 'ntrial', 'Side', and
%                  other velocity/spike/VM fields.
%
%   Outputs:
%       concat_analysis - Struct array where trials are grouped according to
%                         'ntrial', data is concatenated, and 'Side' is assigned.
%
%   Example:
%       concat_analysis = groupAndConcatAnalysis(analysis);

    %% Step 1: Create groups based on ntrial
    groups = {};                % Cell array to store groups
    currentGroup = [];          % Temporary variable to hold current group
    currentBinary = [];         % Temporary variable to hold Side values (binary)
    
    data = [analysis.ntrial]';
    side = [analysis.Side]';
    side = strcmp(side, 'left'); % 1 for left, 0 for right
    
    for i = 1:length(data)
        if data(i) == 1
            if ~isempty(currentGroup)
                groups{end+1, 1} = currentGroup;
                groups{end, 2} = currentBinary;
                currentGroup = [];
                currentBinary = [];
            end
        end
        currentGroup = [currentGroup, i];
        currentBinary = [currentBinary, side(i)];
    end
    
    % Add the last group
    if ~isempty(currentGroup)
        groups{end+1, 1} = currentGroup;
        groups{end, 2} = currentBinary;
    end

    %% Step 2: Define fields to concatenate
    list_of_columns = {'VM'; 'spikesbinned'; 'zveloc_in_mm_binned'; ...
        'zveloc_in_degree_per_s_binned'; 'xveloc_in_mm_binned'; 'VM_medfilt'; ...
        'xveloc_in_mm'; 'zveloc_in_degree_per_s'; 'Side'};

    %% Step 3: Ensure all fields are column vectors
    for k = 1:length(list_of_columns)
        for i = 1:length(analysis)
            if size(analysis(i).(list_of_columns{k}), 1) == 1
                analysis(i).(list_of_columns{k}) = analysis(i).(list_of_columns{k})';
            end
        end
    end

    %% Step 4: Initialize output struct correctly
    concat_analysis = struct();
    for i = 1:length(groups)
        for k = 1:numel(list_of_columns)
            concat_analysis(i).(list_of_columns{k}) = [];  % initialize each field per struct element
        end
    end


    %% Step 5: Concatenate data per group
    for k = 1:length(list_of_columns)
        for i = 1:length(groups)
            rowIndices = groups{i, 1};
            concat_analysis(i).(list_of_columns{k}) = vertcat(analysis(rowIndices).(list_of_columns{k}));
        end
    end

    %% Step 6: Assign Side to each group
    for i = 1:length(groups)
        if sum(groups{i,2}) > 0
            concat_analysis(i).Side = 'left';
        else
            concat_analysis(i).Side = 'right';
        end
    end
end
