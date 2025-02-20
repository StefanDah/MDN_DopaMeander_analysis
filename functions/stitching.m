
% Initialize variables
groups = {}; % Cell array to store groups
currentGroup = []; % Temporary variable to hold the current group
currentBinary = []; % Temporary variable to hold the binary values of the current group

data = [analysis.ntrial]';
side = [analysis.Side]';
side = strcmp(side, 'left'); % 1 for "left", 0 for "right"

for i = 1:length(data)
    if data(i) == 1
        % If we encounter a 1 and currentGroup is not empty, save it
        if ~isempty(currentGroup)
            groups{end+1, 1} = currentGroup; 
            groups{end, 2} = currentBinary;
            currentGroup = []; % Reset currentGroup
            currentBinary = []; % Reset currentBinary
        end
    end
    % Add the current index and binary value to the currentGroup and currentBinary
    currentGroup = [currentGroup, i];
    currentBinary = [currentBinary, side(i)];
end

% Add the last group if it's not empty
if ~isempty(currentGroup)
    groups{end+1, 1} = currentGroup;
    groups{end, 2} = currentBinary;
end

clearvars concat_analysis list_of_columns

list_of_columns = {'VM';'spikesbinned',; 'zveloc_in_mm_binned'; 'zveloc_in_degree_per_s_binned';...
    'xveloc_in_mm_binned'; 'VM_medfilt'; 'xveloc_in_mm'; 'zveloc_in_degree_per_s'; 'Side'};

concat_analysis = struct();

for k = 1 : length(list_of_columns)
    
    
    for i = 1:length(analysis)
        % Check the size of the specific column (e.g., 'VM') for each row
        if size(analysis(i).(list_of_columns{k}), 1) == 1  % If it's a row vector (1 x n)
            % Transform it to a column vector (n x 1)
            analysis(i).(list_of_columns{k}) = analysis(i).(list_of_columns{k})';
        end
    end
    
end
% Loop through column names and initialize each field
for k = 1:numel(list_of_columns)
    concat_analysis.(list_of_columns{k}) = [];  % Assign an empty array
end

for k = 1:length(list_of_columns)
    
    %     concat_analysis.(list_of_columns{k}) = [];
    
    for i = 1:length(groups)
        
        rowIndices = groups{i,1};
    
        
        concat_analysis(i).(list_of_columns{k}) = vertcat(analysis(rowIndices).(list_of_columns{k}));

    end
end

for i = 1 : length(groups)
    if sum(groups{i,2}) > 0
        concat_analysis(i).Side = 'left';
    else
        concat_analysis(i).Side = 'right';
    end
end
