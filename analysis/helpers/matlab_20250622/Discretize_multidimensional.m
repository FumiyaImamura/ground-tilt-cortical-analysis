function [JointDescretizedData, F] = Discretize_multidimensional(Data_anyDimensional, edges_cellVector, dim2jointlyDiscretize)
% DISCRETIZE_MULTIDIMENSIONAL Jointly discretizes a multidimensional coordinate.
%
%   [voxelIdx, F] = Discretize_multidimensional(Data_anyDimensional, edges_cellVector, dim2jointlyDiscretize)
%
%   Input:
%     Data_anyDimensional   - A multidimensional array where the coordinate dimension is given by 
%              the scalar 'dim2jointlyDiscretize'. That is, along dimension 'dim2jointlyDiscretize' the coordinate values
%              for each datapoint are stored; the other dimensions index datapoints.
%
%     edges_cellVector  - A cell vector of length equal to size(Data_anyDimensional,dim2jointlyDiscretize). Each cell holds a vector
%              of bin edges_cellVector to use for the corresponding coordinate.
%
%     dim2jointlyDiscretize    - The dimension of Data_anyDimensional that corresponds to the coordinate.
%
%   Output:
%     JointDescretizedData - An array (with the coordinate dimension removed, i.e. replaced by 1)
%							containing for each datapoint a single linear voxel index. Its shape 
%							is the same as Data_anyDimensional except that size(voxelIdx,dim2jointlyDiscretize)==1.
%
%     F        - A structure with fields:
%                  F.discretize(newdata) returns the voxel indices for newdata (which 
%                                         must have the coordinate dimension in position 'dim2jointlyDiscretize')
%                  F.undiscretize(discData, bincenter) returns for each voxel index the 
%                                         voxel center coordinates. The output has the 
%                                         coordinate dimension inserted at position 'dim2jointlyDiscretize'.
%                  F.edges  (same as input edges_cellVector)
%                  F.dim    (same as input dim2jointlyDiscretize)
%                  F.bincount  = total number of voxels (product of bins per coordinate)
%
%   The voxel centers are computed as follows:
%     - If a nonempty cell in bincenter is provided for a coordinate, its value for a
%       given voxel is taken from that cell.
%     - Otherwise, the voxel center is computed as the midpoint between the lower and
%       upper bin edges_cellVector.
%
%   Example:
%     % Suppose Data_anyDimensional is a 3-by-100 array (each column a 3D coordinate) with coordinate
%     % dimension 1.
%     edges_cellVector = {linspace(0,1,11), linspace(-1,1,21), linspace(100,200,11)};
%     [voxelIdx, F] = Discretize_multidimensional(Data_anyDimensional, edges_cellVector, 1);
%
%     % Discretize new data:
%     newVoxelIdx = F.discretize(newData);
%
%     % Get the voxel center coordinates; here, bincenter is optional.
%     centers = F.undiscretize(newVoxelIdx, bincenter);
%

Data=Data_anyDimensional;
edges=edges_cellVector;
dim=dim2jointlyDiscretize;

    %% --- Input validation ---
    if ~iscell(edges)
        % error('edges must be provided as a cell vector.');
		if isvector(edges)
			edges={edges};	% if it was a numeric vector, just make it into a single cell
		end
    end
    nCoord = size(Data, dim);
    if length(edges) ~= nCoord
        error('Length of edges (%d) must equal size(Data,dim) (%d).', length(edges), nCoord);
    end

    %% --- Permute Data so that the coordinate dimension is first ---
    % Build a permutation so that the coordinate (dim) becomes the first dimension.
    numDims = ndims(Data);
    permOrder = [dim, 1:dim-1, dim+1:numDims];
    DataPerm = permute(Data, permOrder);
    
    %% --- Flatten the datapoint dimensions ---
    sPerm = size(DataPerm);  % sPerm(1)=nCoord, sPerm(2:end)= datapoint dims.
    nPoints = prod(sPerm(2:end));
    DataReshaped = reshape(DataPerm, sPerm(1), nPoints);
    
    %% --- Discretize each coordinate using the provided edges ---
    binIndices = nan(sPerm(1), nPoints);
    nbins = zeros(sPerm(1), 1);  % number of bins per coordinate.
    for i = 1:sPerm(1)
        binIndices(i, :) = discretize(DataReshaped(i, :), edges{i});
        nbins(i) = numel(edges{i}) - 1;
    end
    
    %% --- Combine the multidimensional bin indices into one linear index ---
    % The mapping is based on:
    %   linear_index = 1 + sum( (subscript_i - 1) * multiplier(i) )
    multiplier = ones(sPerm(1), 1);
    for i = 2:sPerm(1)
        multiplier(i) = multiplier(i-1) * nbins(i-1);
    end
    voxelLinear = 1 + sum((binIndices - 1) .* repmat(multiplier, 1, nPoints), 1);
    
    %% --- Reshape the voxel indices to the original datapoint configuration ---
    % We want the output to have the same shape as Data except that the coordinate 
    % dimension is replaced by a singleton.
    origSize = size(Data);
    if numel(origSize) < dim
        origSize(end+1:dim) = 1;
    end
    outSize = origSize;
    outSize(dim) = 1;
    voxelIdx = reshape(voxelLinear, outSize);
    
	JointDescretizedData=voxelIdx;

    %% --- Define function handles in output structure F ---
    F.discretize   = @discretizeNewData;
    F.undiscretize = @undiscretize;
    F.edges        = edges;
    F.dim          = dim;
    F.bincount     = prod(nbins, "all");
    
    %% --- Nested function: discretizeNewData ---
    function idx = discretizeNewData(newdata)
        % DISCRETIZENEWDATA Discretizes new data in the same manner.
        %
        %   idx = F.discretize(newdata)
        %
        %   newdata must have the coordinate dimension in position 'dim'.
        
        if size(newdata, dim) ~= nCoord
            error('newdata must have size(newdata,dim) = %d (matching the number of edges).', nCoord);
        end
        
        % Permute so that coordinate dimension comes first.
        newdataPerm = permute(newdata, permOrder);
        s_new = size(newdataPerm);  % s_new(1)=nCoord; s_new(2:end)= datapoint dims.
        nPoints_new = prod(s_new(2:end));
        newdataReshaped = reshape(newdataPerm, s_new(1), nPoints_new);
        
        % Discretize each coordinate.
        binIdx = nan(s_new(1), nPoints_new);
        for j = 1:s_new(1)
            binIdx(j, :) = discretize(newdataReshaped(j, :), edges{j});
        end
        
        % Combine the bin indices into a linear index.
        idxLinear = 1 + sum((binIdx - 1) .* repmat(multiplier, 1, nPoints_new), 1);
        
        % Reshape so that the output has the same shape as newdata except that
        % the coordinate dimension is replaced by a singleton.
        origSize_new = size(newdata);
        if numel(origSize_new) < dim
            origSize_new(end+1:dim) = 1;
        end
        outSize_new = origSize_new;
        outSize_new(dim) = 1;
        idx = reshape(idxLinear, outSize_new);
    end

    %% --- Nested function: undiscretize ---
    function centers = undiscretize(discData, bincenter)
        % UNDISCRETIZE Maps voxel indices back to voxel center coordinates.
        %
        %   centers = F.undiscretize(discData, bincenter)
        %
        %   discData is an array of voxel indices (for example, produced by F.discretize).
        %   
		%   bincenter is an optional cell array (of length nCoord) where each cell holds a vector
        %   of bin centers for that coordinate. For any coordinate, if bincenter{k} is empty (or if
        %   bincenter itself is empty), the voxel center is computed as the midpoint of the bin edges.
        %
        %   The output 'centers' is an array with the coordinate dimension inserted at position 'dim'.
        %   That is, if discData has size S (which corresponds to the non-coordinate dimensions),
        %   centers will have size:
        %         [S(1:dim-1), nCoord, S(dim:end)].
        
        % Decide whether to use bincenter.
        if nargin < 2 || isempty(bincenter)
            % useBincenter = false;
			%%%%%%%%%
			useBincenter=true;
			bincenter=cellfun(@(a)make_bins_from_edges(a),edges,UniformOutput=false);
        else
            useBincenter = true;
            if ~iscell(bincenter)
                error('bincenter must be a cell array.');
            end
            if length(bincenter) ~= nCoord
                error('Length of bincenter (%d) must match the number of coordinates (%d).', length(bincenter), nCoord);
            end
        end
        
        % Flatten discData.
        L = discData(:);
        nL = numel(L);
        % Initialize an array to hold the bin (subscript) indices for each coordinate.
        subs = nan(nCoord, nL);
        
        % Since we assume no NaN values occur, process all indices.
        temp = L - 1;
        for k = 1:nCoord
            subs(k, :) = mod(temp, nbins(k)) + 1;
            temp = floor(temp / nbins(k));
        end
        
        % Compute the voxel centers for each coordinate.
        centersVals = nan(nCoord, nL);
        for k = 1:nCoord
            idxs = subs(k, :);
            if useBincenter && ~isempty(bincenter{k})
                centersVals(k, :) = bincenter{k}(idxs);
            else
                % Compute the midpoint of the bin edges.
                centersVals(k, :) = (edges{k}(idxs) + edges{k}(idxs+1)) / 2;
            end
        end
        
        % Reshape and insert the coordinate dimension at position 'dim'.
        % We want the output shape to be the same as discData except that the
        % coordinate dimension has size nCoord.
        % Start by reshaping each coordinate's values to the shape of discData.
        centers = [];
        for k = 1:nCoord
            % reshape to the same shape as discData (which is the same as voxelIdx output)
            curCoord = reshape(centersVals(k, :), size(discData));
            % Concatenate along dimension 'dim'.
            centers = cat(dim, centers, curCoord);
        end
    end

end
