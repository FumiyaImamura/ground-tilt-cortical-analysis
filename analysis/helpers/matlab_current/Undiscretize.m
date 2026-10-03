function undiscretized = Undiscretize(discretized, edges)
% undiscretize_array Converts a discretized array back to its undiscretized values.
%
%   undiscretized = undiscretize_array(discretized, edges)
%
%   This function takes a multi-dimensional discretized array (where each element is
%   assumed to be a bin index) and a vector of bin edges. For each bin, it calculates
%   the representative value as the mean of the two edges that form the bin.
%
%   If either the lower or upper edge is infinite (i.e. -Inf or Inf), these are replaced
%   by finite values based on the mean interval calculated from the finite edges.
%
%   In addition, if any element in the discretized array is NaN or not a valid bin index 
%   (i.e. it is not a positive integer between 1 and (numel(edges)-1)), the function 
%   outputs a warning. The corresponding elements in the output will be set to NaN.
%
%   Inputs:
%       discretized - A multi-dimensional array containing bin indices.
%       edges       - A vector containing the edges used for discretization. This vector
%                     may include -Inf or Inf.
%
%   Output:
%       undiscretized - An array of the same size as "discretized" where each valid bin
%                       index has been replaced by its representative (mean) value.
%
%   Example:
%       edges = [-Inf, 0, 1, Inf];
%       dArray = [1, 2, 1; 2, 1, 2];
%       undisArray = undiscretize_array(dArray, edges);
%

    % Calculate bin centers based on the given edges.
    % First, isolate the finite edges.
    finiteEdges = edges(~isinf(edges));
    if numel(finiteEdges) < 2
        error('Not enough finite edges to compute a mean interval.');
    end
    
    % Compute the mean interval of the finite edges.
    meanInterval = mean(diff(finiteEdges));

    % Replace -Inf and Inf with a finite value.
    newEdges = edges;
    if isinf(newEdges(1))
        newEdges(1) = finiteEdges(1) - meanInterval;
    end
    if isinf(newEdges(end))
        newEdges(end) = finiteEdges(end) + meanInterval;
    end

    % Calculate the bin centers as the mean of consecutive edges.
    binCenters = (newEdges(1:end-1) + newEdges(2:end)) / 2;

    % Prepare the output array. Start with NaN to mark invalid entries.
    undiscretized = nan(size(discretized));

    % Check for any invalid bin indices.
    % A valid bin index is a positive integer from 1 to numel(binCenters).
    invalidMask = isnan(discretized) | ...
                  (discretized < 1) | ...
                  (discretized > numel(binCenters)) | ...
                  (mod(discretized, 1) ~= 0);
    if any(invalidMask(:))
        warning('The discretized array contains NaN or invalid bin indices. These elements will be set to NaN in the output.');
    end

    % For each valid bin index, assign the corresponding bin center value.
    validMask = ~invalidMask;
    for i = 1:numel(binCenters)
        % Only update the indices that exactly equal the bin index i.
        idx = validMask & (discretized == i);
        undiscretized(idx) = binCenters(i);
    end
end
