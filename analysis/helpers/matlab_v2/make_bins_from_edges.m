function binCenters = make_bins_from_edges(edges)
% make_bins_from_edges Returns the representative bin values from a vector of edges.
%
%   binCenters = make_bins_from_edges(edges)
%
%   This function takes a vector of bin edges, which may include -Inf or Inf, and
%   computes the bin centers. Each bin center is calculated as the average of two
%   consecutive edge values. If the first or last edge is infinite, it is replaced by
%   a finite value determined by the mean interval computed from the finite edges.
%
%   Inputs:
%       edges - A vector containing the bin edges (for example, [-Inf, 0, 1, Inf]).
%
%   Output:
%       binCenters - A vector containing the representative value (center) for each bin.
%
%   Example:
%       edges = [-Inf, 0, 1, Inf];
%       centers = make_bins_from_edges(edges);
%       % centers will be a vector containing the center value for each bin.
%

    % Identify the finite edges for calculating the mean interval.
    finiteEdges = edges(~isinf(edges));
    if numel(finiteEdges) < 2
        error('Not enough finite edges to compute a mean interval.');
    end

    % Compute the mean interval using the differences between finite edges.
    meanInterval = mean(diff(finiteEdges));

    % Create a modified version of the edges vector.
    newEdges = edges;
    if isinf(newEdges(1))
        % Replace -Inf with the first finite edge minus the mean interval.
        newEdges(1) = finiteEdges(1) - meanInterval;
    end
    if isinf(newEdges(end))
        % Replace Inf with the last finite edge plus the mean interval.
        newEdges(end) = finiteEdges(end) + meanInterval;
    end

    % Compute bin centers: the representative value is the mean of consecutive edges.
    binCenters = (newEdges(1:end-1) + newEdges(2:end)) / 2;
end
