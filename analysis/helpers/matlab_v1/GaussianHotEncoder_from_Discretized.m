function G=GaussianHotEncoder_from_Discretized(Data_numeric,encoding_dim,parameter_dim_mustBeSize1,sigma,bincount)
% GaussianHotEncoder_from_Discretized applies a Gaussian smoothing to discretized data.
% It converts the discretized values along the specified parameter dimension 
% into Gaussian curves along the encoding dimension (typically time). The smoothing width is 
% controlled by sigma and the number of discrete bins by bincount.
%
% Input:
%   Data_numeric         - Numeric array with discretized data
%   encoding_dim - Dimension along which the Gaussian is applied
%   parameter_dim_mustBeSize1- Dimension that encoded parameters are going to align in (must be size 1)
%   sigma        - Standard deviation of the Gaussian curve
%   bincount     - Number of bins for the discretization
%
% Output:
%   Data_disc_numeric    - The Gaussian-smoothed data array

arguments
    Data_numeric            % numeric array (discretized). size(Data_numeric,parameter_dim_mustBeSize1) must be 1 since this is the dimension that gaus hot encoding parameters are going to be
    encoding_dim
    parameter_dim_mustBeSize1    
    sigma (1,:)
    bincount
end

% Build Gaussian kernel (±3σ)
w = ceil(3*sigma);
g = exp(-0.5*((-w:w)/sigma).^2)';
g = g./sum(g);

G=HotEncoder(Data_numeric,encoding_dim,parameter_dim_mustBeSize1,g,bincount);

if false
	%%
	% --- Example script for GaussianHotEncoder ---

	% Parameters
	T       = 100;    % # time points
	nBins   = 5;      % # discrete bins
	nTrials = 2;      % # “trials” or repetitions
	sigma   = 0.5;      % width of the Gaussian

	% 1) Create synthetic discretized data:
	%    Data has size [time × 1 × trials], with values in 1…nBins
	Data = randi(nBins, [T, 1, nTrials]);

	% 2) Apply the encoder
	%    encoding_dim = 1 (time), param_dim = 2 (singleton)
	G = GaussianHotEncoder_from_Discretized(Data, 1, 2, sigma, nBins);

	% 3) Prepare time vector and colors
	time   = (1:T)';
	colors = lines(nBins);

	% 4) Plot the Gaussian curves for trial #1
	figure; hold on
	for b = 1:nBins
		plot(time, squeeze(G(:,b,1)), 'Color', colors(b,:), 'LineWidth', 1.5);
	end

	% 5) Overlay original one‑hot impulses with matching colors
	for b = 1:nBins
		idx = (Data(:,1,1) == b);
		scatter(time(idx), zeros(sum(idx),1) + 0.02, ...
			40, colors(b,:), 'filled', 'MarkerEdgeColor', 'none');
	end

	xlabel('Time');
	ylabel('Amplitude');
	legend(arrayfun(@(x) sprintf('Bin %d',x), 1:nBins, 'UniformOutput', false), ...
		'Location','best');
	title('Gaussian Hot Encoding (trial 1)');
	hold off
end
end