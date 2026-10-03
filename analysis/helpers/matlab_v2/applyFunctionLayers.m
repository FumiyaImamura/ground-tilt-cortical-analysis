function output = applyFunctionLayers(inputVar, funcHandle,options)
arguments
	inputVar
	funcHandle
	options.AcceptStrings (1,1) logical =true
end
    % applyFunctionLayers applies a function to all numeric elements
    % within a nested structure of cell arrays, structs, or numeric arrays.
    % * The deepest elements that are not numeric are skipped!
    % Inputs:
    %   - inputVar: A numeric array, cell array, or struct (possibly nested)
    %   - funcHandle: A function handle to be applied to numeric elements
    %
    % Output:
    %   - output: The modified structure with the function applied


    % Check if the input is numeric
	if isnumeric(inputVar) || (options.AcceptStrings & isstring(inputVar))
        % Apply the function handle directly
        output = funcHandle(inputVar);

    % Check if the input is a cell array
    elseif iscell(inputVar)
		options=namedargs2cell(options);
        % Initialize the output as a cell array of the same size
        output = cell(size(inputVar));

        % Iterate over each element in the cell array
        for i = 1:numel(inputVar)
            % Recursively call the function for each element
            output{i} = applyFunctionLayers(inputVar{i}, funcHandle,options{:});
        end

    % Check if the input is a struct
    elseif isstruct(inputVar)
		options=namedargs2cell(options);
        if numel(inputVar) > 1
            % If the struct is an array (non-scalar), process each element individually.
            output = inputVar;  % Preallocate the output with the same size.
            for i = 1:numel(inputVar)
                output(i) = applyFunctionLayers(inputVar(i), funcHandle,options{:});
            end
        else
            % For a scalar struct, process field by field.
            output = struct();
            flds = fieldnames(inputVar);
            for i = 1:numel(flds)
                field = flds{i};
                % Recursively call the function for each field value.
                output.(field) = applyFunctionLayers(inputVar.(field), funcHandle,options{:});
            end
        end

    else
        % If the input is neither numeric, cell array, nor struct, skip it.
        output = inputVar;
    end
end
