function CellWithStruct = StructWithCell2CellWithStruct_multiLayers(StructWithCell)
% STRUCTWITHCELL2CELLWITHSTRUCT Recursively converts a structure whose 
% fields contain cell arrays into a cell array of structures.
%
%   CellWithStruct = StructWithCell2CellWithStruct(StructWithCell)
%
% The function searches through the input structure, and whenever it
% finds a cell array (even if nested inside another structure) it
% expands that layer. Scalar fields are replicated to match the expanded
% dimensions. In the example:
%
%   StrucWithCell.A.B = {[1,2], 2, 3};
%
% the function produces:
%
%   CellWithStruct{1}.A.B = [1,2];
%   CellWithStruct{2}.A.B = 2;
%   CellWithStruct{3}.A.B = 3;
%
% If no cell array is encountered at the top level, the final result is
% wrapped into a cell array.

    % Use the recursive helper function
    result = expandStruct(StructWithCell);
    
    % Ensure that the output is a cell array.
    if ~iscell(result)
        CellWithStruct = {result};
    else
        CellWithStruct = result;
	end


	if false
		%% Example code to assess the code accuracy
		A=[];
		A.A={1,2,3,4};
		A.B.a={2,3,4,5};
		A.B.b={3,4,5,6};
		A.C.a.a={4,5,6,7};
		newA=StructWithCell2CellWithStruct_multiLayers(A);
	end
end

function out = expandStruct(s)
% EXPANDSTRUCT Recursively processes a structure.
%
% If any field of the structure is a cell array, then the structure is
% expanded into a cell array of structures, one per element of the cell.
% Scalar (non-cell) fields are replicated for each expanded structure.

    if ~isstruct(s)
        % If input is not a structure, no further expansion is possible.
        out = s;
        return;
    end

    % Get all field names in the structure.
    flds = fieldnames(s);
    
    % Process each field recursively.
    for i = 1:length(flds)
        f = flds{i};
        s.(f) = expandField(s.(f));
    end
    
    % Look for fields that are cell arrays.
    cellFields = {};
    len = [];
    for i = 1:length(flds)
        f = flds{i};
        if iscell(s.(f))
            cellFields{end+1} = f;  %#ok<AGROW>
            len(end+1) = numel(s.(f)); %#ok<AGROW>
        end
    end
    
    if ~isempty(cellFields)
        % If one or more fields are cell arrays, assume they have the same length.
        if any(len ~= len(1))
            error('Inconsistent cell array lengths among fields.');
        end
        N = len(1);
        out = cell(1, N);
        % For each element in the cell array(s), create a copy of the structure.
        for j = 1:N
            temp = s; % start from the full structure copy
            for k = 1:length(cellFields)
                fieldName = cellFields{k};
                % Replace cell fields with their j-th element.
                temp.(fieldName) = s.(fieldName){j};
            end
            % In case the replacement inserted another cell array (nested),
            % process the structure again.
            out{j} = expandStruct(temp);
        end
    else
        % If no cell array is found, return the structure unchanged.
        out = s;
    end
end

function out = expandField(v)
% EXPANDFIELD Processes an individual field value.
%
% If the value is a structure, it is processed via expandStruct.
% If it is a cell array, each element is processed recursively.
% Otherwise, the value is left unchanged.

    if isstruct(v)
        out = expandStruct(v);
    elseif iscell(v)
        % Process each element in the cell array.
        out = cellfun(@expandField, v, 'UniformOutput', false);
    else
        out = v;
    end
end

