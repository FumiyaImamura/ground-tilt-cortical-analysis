function hasNaN = assessNaN(inputVar)
    % assessNaN checks if there are any NaN values within a nested structure
    % of numeric arrays, cell arrays, or structs.
    %
    % Input:
    %   - inputVar: A numeric array, cell array, or struct (possibly nested)
    %
    % Output:
    %   - hasNaN: A boolean that is true if any numeric element contains NaN,
    %             false otherwise.
    
    % Use applyFunctionLayers with a function handle that checks for NaN in a numeric array.
    % The anonymous function @(x) any(isnan(x(:))) returns true if any element of x is NaN.
    nanFlags = applyFunctionLayers(inputVar, @(x) any(isnan(x(:))));
    
    % Now, nanFlags is a nested structure (or a scalar) containing booleans.
    % We need to aggregate these booleans into one overall result.
    hasNaN = aggregateNaNFlags(nanFlags);
end

function flag = aggregateNaNFlags(data)
    % aggregateNaNFlags recursively checks a nested structure of booleans
    % (or mixed types) and returns true if any true is found.
    if islogical(data)
        % Base case: if data is logical, return its value.
        flag = data;
    elseif iscell(data)
        % If data is a cell array, iterate over each element.
        flag = false;
        for i = 1:numel(data)
            flag = flag || aggregateNaNFlags(data{i});
            if flag
                break;  % short-circuit if true found
            end
        end
    elseif isstruct(data)
        % If data is a struct, iterate over each field.
        flag = false;
        fields = fieldnames(data);
        for i = 1:numel(fields)
            flag = flag || aggregateNaNFlags(data.(fields{i}));
            if flag
                break;  % short-circuit if true found
            end
        end
    else
        % For any other type (shouldn't happen here), return false.
        flag = false;
    end
end
