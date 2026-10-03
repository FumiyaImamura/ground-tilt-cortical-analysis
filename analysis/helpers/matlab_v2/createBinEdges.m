function [edges_numeric_vec,bincenter_numeric_vec] = createBinEdges(data_numeric, bin_count, method, percentile_threshold_for_Uniform_closeTo0)
% createBinEdges generates bin edges for data binning.
%
% This function takes in a numeric array and creates edges to be used for binning the
% data into a specified number of bins. There are two methods to compute these edges:
%
%   1. 'Uniform'
%      In the Uniform method, the user must provide a percentile_threshold_for_Uniform (a scalar between
%      0 and 100). The central range is defined by the percentile_threshold_for_Uniform and (100 - percentile_threshold_for_Uniform)
%      percentiles of the data. The inner bin edges between these percentiles are uniformly spaced.
%      The full set of edges is constructed by prepending -Inf and appending Inf. This means that all
%      data below the lower percentile will fall into the first bin, and all data above the upper
%      percentile will fall into the last bin.
%
%   2. 'Proportional'
%      In the Proportional method, the edges are computed so that each bin contains roughly the same
%      number of data points. The bin edges are determined by computing the percentiles that partition
%      the data into equal proportions. In this case, the percentile_threshold_for_Uniform argument is ignored.
%
% Syntax:
%   edges = createBinEdges(data, bin_count, method, percentile_threshold_for_Uniform)
%
% Inputs:
%   data_numeric               - Numeric array of any dimension.
%   bin_count          - A positive integer indicating the number of bins.
%   method             - A string specifying the method; must be either 'Uniform' or 'Proportional'.
%   percentile_threshold_for_Uniform - (Required only for 'Uniform') A scalar value between 0 and 100 that specifies
%                          the percentiles used to define the lower and upper bounds of the central range.
%
% Outputs:
%   edges - A vector of bin edges:
%           * For 'Uniform': edges = [-Inf, inner_edge_1, inner_edge_2, ..., inner_edge_n, Inf]
%           * For 'Proportional': edges are the data values at the corresponding percentiles,
%             spanning from the 0th to the 100th percentile.
%
% Examples:
%   % Example using 'Uniform' method with 4 bins and a 5% threshold:
%   data = randn(1000,1);
%   edges = createBinEdges(data, 4, 'Uniform', 5);
%
%   % Example using 'Proportional' method with 4 bins:
%   edges = createBinEdges(data, 4, 'Proportional');
%
% Notes:
%   - For 'Uniform' method, the range defined by the specified percentiles is used to determine
%     inner edges; data outside this range is grouped into the first or last bin.
%   - The function uses MATLAB's "arguments" block for input validation.
%


    % Validate inputs using the arguments block.
    arguments
        data_numeric {mustBeNumeric}
        bin_count (1,1) {mustBeInteger, mustBePositive}
        method (1,:) char {mustBeMember(method, {'Uniform', 'Proportional'})}="Uniform"
        percentile_threshold_for_Uniform_closeTo0 double = NaN  % Default value; will be validated if needed.
	end
	data=data_numeric;
	percentile_threshold_for_Uniform=percentile_threshold_for_Uniform_closeTo0;

    % For the 'Uniform' method, ensure that a valid percentile_threshold_for_Uniform is provided.
    if strcmp(method, 'Uniform')
        if isnan(percentile_threshold_for_Uniform)
            error('For "Uniform" method, a valid percentile_threshold_for_Uniform (between 0 and 100) must be provided.');
        else
            % Confirm that the percentile_threshold_for_Uniform is within the acceptable range.
            validateattributes(percentile_threshold_for_Uniform, {'numeric'}, {'scalar', '>=', 0, '<=', 100});
        end
    end

    % Flatten the input data to a column vector for easier percentile computations.
    data = data(:);

    % Compute bin edges based on the chosen method.
    if strcmp(method, 'Uniform')
        % -----------------------------
        % Uniform Method Explanation:
        %   1. Compute the lower bound as the "percentile_threshold_for_Uniform"-th percentile.
        %   2. Compute the upper bound as the (100 - percentile_threshold_for_Uniform)-th percentile.
        %   3. If more than one bin is requested, generate uniformly spaced inner edges between
        %      the lower and upper bounds.
        %   4. Prepend -Inf and append Inf so that data outside the central range is included
        %      in the first and last bins.
        % -----------------------------
        
        % Determine the lower and upper bounds from the data.
        lower_bound = prctile(data, percentile_threshold_for_Uniform);
        upper_bound = prctile(data, 100 - percentile_threshold_for_Uniform);
        
        % Special case: if only one bin is requested.
        if bin_count == 1
            edges_numeric_vec = [-Inf, Inf];
        else
            % Generate uniformly spaced inner edges.
			%%%%%%%%%%%%%%%% update 2025/04/18	
			% Previously, edges were selected so that edge bins contained
			% only percentile_threshold_for_Uniform % of data. 
			% Now, it is percentile_threshold_for_Uniform % + the width of
			% bins
			inner_edges = linspace(lower_bound, upper_bound, bin_count + 1);
			inner_edges= inner_edges(2:end-1);
			%%%%%%%%%%%%%%%
            % inner_edges = linspace(lower_bound, upper_bound, bin_count - 1);
			%%%%%%%%%%%%%%%%
            % Construct the complete set of edges.
            edges_numeric_vec = [-Inf, inner_edges, Inf];
        end

    elseif strcmp(method, 'Proportional')
        % -----------------------------
        % Proportional Method Explanation:
        %   1. The function computes equally spaced percentiles that partition the data into bins
        %      with roughly the same number of data points.
        %   2. For bin_count bins, (bin_count+1) percentiles are computed, from 0 to 100.
        % -----------------------------
        
        % Create percentiles for the bin edges.
        percents = linspace(0, 100, bin_count + 1);
        edges_numeric_vec = prctile(data, percents);
	end

	bincenter_numeric_vec=diff(edges_numeric_vec,1)/2+edges_numeric_vec(1:end-1);

end
