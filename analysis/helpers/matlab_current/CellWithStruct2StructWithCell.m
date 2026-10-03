function outStruct = CellWithStruct2StructWithCell(cellArray)
% CellWithStruct2StructWithCell Converts a cell array of structs into a struct of cells.
%
%   outStruct = CellWithStruct2StructWithCell(cellArray)
%
%   This function takes a cell array where each element is a struct (with possibly
%   different fields) and produces a new struct. The new struct has fields corresponding
%   to the union of all fields from all the structs in the input cell array. Each field
%   in the output struct is a cell array of the same size as the input cell array.
%   If an element in the input cell array does not contain a particular field, the corresponding
%   cell in the output struct's field will be empty.
%
%   Example:
%       % Define a cell array of structs with different fields:
%       C = {struct('a', 1, 'b', 2), struct('b', 3, 'c', 4)};
%
%       % Convert to a struct with cell arrays:
%       S = CellWithStruct2StructWithCell(C);
%
%       % Now, S.a = {1, []}, S.b = {2, 3}, S.c = {[], 4}.
%

    % Validate input is a cell array
    if ~iscell(cellArray)
        error('Input must be a cell array.');
    end

    % Determine the union of all field names from each struct in the cell array.
    allFields = {};
    for idx = 1:numel(cellArray)
        currStruct = cellArray{idx};
        if ~isempty(currStruct)
            if ~isstruct(currStruct)
                error('Each element of the cell array must be a struct.');
            end
            % Get the fieldnames from the current struct and merge with allFields.
            currFields = fieldnames(currStruct);
            allFields = union(allFields, currFields);
        end
    end

    % Initialize the output struct. For each field in the union, create a cell array
    % of the same size as the input cell array.
    outStruct = struct();
    cellSize = size(cellArray);
    for i = 1:length(allFields)
        outStruct.(allFields{i}) = cell(cellSize);
    end

    % Populate the output struct: for each cell element in the input,
    % if the struct has a given field, assign its value to the corresponding cell in the output.
    for idx = 1:numel(cellArray)
        currStruct = cellArray{idx};
        if isempty(currStruct)
            continue;
        end
        currFields = fieldnames(currStruct);
        for j = 1:length(currFields)
            fieldName = currFields{j};
            outStruct.(fieldName){idx} = currStruct.(fieldName);
        end
    end

end
