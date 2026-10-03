function Data_disc = GaussianHotEncoder(Data, encoding_dim, parameter_dim,sigma, edges)
% GaussianHotEncoder - Apply Gaussian hot encoding along the time dimension for each parameter and trial.
%
% This function discretizes each parameter (located along the specified parameter dimension of Data)
% by first mapping the data to discrete bins using the provided edges and then replacing the resulting
% one-hot encoding with a Gaussian curve along the encoding (typically time) dimension.
%
% The function accepts two types of 'edges' inputs:
%    1) **Cell array input**: Each cell contains a vector of bin edges for the corresponding parameter.
%       In this case, the output Data_disc is returned as a cell array, with each cell containing the
%       Gaussian-hot encoded data for that parameter. If a cell in the 'edges' array is empty, that
%       parameter is left unprocessed.
%
%    2) **Non-cell input** (e.g., a numeric vector): The numeric vector is used for discretizing a
%       single parameter. In this case, Data_disc is returned as a numeric array.
%
% For parameters that are discretized, the second dimension of the output is expanded to reflect the
% number of bins.
%
% Syntax:
%    Data_disc = GaussianHotEncoder(Data, encoding_dim, parameter_dim, sigma, edges)
%
% Inputs:
%    Data          - A multi-dimensional array containing the data to be processed.
%    encoding_dim  - The dimension index corresponding to the time (or encoding) points along which
%                    the Gaussian encoding is applied.
%    parameter_dim - The dimension index corresponding to the parameters that need to be Gaussian-hot encoded.
%    sigma         - The standard deviation of the Gaussian curve used to "soften" the one-hot encoding.
%    edges         - Either a cell array where each cell contains a vector of bin edges for the corresponding
%                    parameter, or a numeric vector for discretizing a single parameter.
%
% Outputs:
%    Data_disc     - If 'edges' is a cell array, a cell array containing the Gaussian-hot encoded outputs
%                    for each parameter (with the second dimension expanded to reflect the number of bins).
%                    If 'edges' is a numeric vector, a numeric array containing the Gaussian-hot encoded data.
%
% Example:
%    % Example 1: Multiple parameters using a cell array for edges
%    edges = {linspace(0,1,11), [], linspace(-2,2,9)};
%    Data_disc = GaussianHotEncoder(Data, 1, 2, sigma, edges);
%
%    % Example 2: Single parameter using a numeric vector for edges
%    edges = linspace(0,1,11);
%    Data_disc = GaussianHotEncoder(Data, 1, 2, sigma, edges);
%

%---------------------------------------------------------------------------
arguments
    Data
	encoding_dim
	parameter_dim
    sigma
    edges
end
%% move the encoding dimension to 1st and parameter dimension to 2nd dimension
v=1:max([ndims(Data),parameter_dim,encoding_dim]);
firstdims=[encoding_dim,parameter_dim];
laterdims=setdiff(v,firstdims);
pvec=[firstdims,laterdims];

orgData=Data;
Data=permute(Data,pvec);


%% concatenate all the higher dimensions to trial dimension
nDims=ndims(Data);
orgsize=size(Data);
if nDims>3
	Data=reshape(Data,size(Data,1),size(Data,2),[]);
end

if iscell(edges)
	if size(Data,2)~=length(edges)
		warning('The number of parameters in Data does not match length of edge. Proceeding anyways...')
	end
else
	if size(Data,2)~=1
		error("The number of parameters and the number of edges do not match");
	end
end



%% generate gaussian
if iscell(edges)
	for p=1:length(edges)
		if ~isempty(edges{p})
			current_edges = edges{p};
			bin_indices = discretize(Data(:, p, :), current_edges);
			if any(isnan(bin_indices(:)))
				warning(append('There are NaN values for edge # ',num2str(p)));
			end
			% find elements that match
			bin_log=(bin_indices-[1:length(current_edges)-1])==0;	% time x parameters (one hot encoded) x trials
			bin_coord=bin_log.*[1:size(bin_log,1)]';
			bin_coord(bin_coord==0)=NaN;
			bin_coord=permute(bin_coord,[4,2,1,3]);					% 1 x parameters (one hot encoded) x time x trials
			indices=[1:size(bin_log,1)]'+zeros(size(bin_log));		% time x parameters (one hot encoded) x trials
			indices=InsertSingleDim(indices,3);						% time x parameters (one hot encoded) x 1 x trials
			gaus_separated = exp(-0.5 * ((indices-bin_coord) / sigma).^2);	% time x parameters (one hot encoded) x time x trials
			gaus_sum=sum(gaus_separated,3,"omitnan");
			gaus_sum=squeeze(gaus_sum);								% time x parameters (one hot encoded) x trials
			if nDims>3
				orgsize(2)=size(gaus_sum,2);
				gaus_sum=reshape(gaus_sum,orgsize);
			end
			Data_disc{p}=gaus_sum;
		else
			Data_disc{p}=[];
		end
		% permute back to the original size
		Data_disc{p}=ipermute(Data_disc{p},pvec);
	end
else
	current_edges = edges;
	bin_indices = discretize(Data, current_edges);
	if any(isnan(bin_indices(:)))
		warning(append('There are NaN values for edge # ',num2str(p)));
	end
	% find elements that match
	bin_log=(bin_indices-[1:length(current_edges)-1])==0;	% time x parameters (one hot encoded) x trials
	bin_coord=bin_log.*[1:size(bin_log,1)]';
	bin_coord(bin_coord==0)=NaN;
	bin_coord=permute(bin_coord,[4,2,1,3]);					% 1 x parameters (one hot encoded) x time x trials
	indices=[1:size(bin_log,1)]'+zeros(size(bin_log));		% time x parameters (one hot encoded) x trials
	indices=InsertSingleDim(indices,3);						% time x parameters (one hot encoded) x 1 x trials
	gaus_separated = exp(-0.5 * ((indices-bin_coord) / sigma).^2);	% time x parameters (one hot encoded) x time x trials
	gaus_sum=sum(gaus_separated,3,"omitnan");
	gaus_sum=squeeze(gaus_sum);								% time x parameters (one hot encoded) x trials
	if nDims>3
		orgsize(2)=size(gaus_sum,2);
		gaus_sum=reshape(gaus_sum,orgsize);
	end
	Data_disc=gaus_sum;
	% permute back to the original size
	Data_disc=ipermute(Data_disc,pvec);
end
end



%%
%%%%%%%%%% RECOMMENDED :
% plotGaussianHotEncodedData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function justPlotThis(Data_disc,edges)
%%

figure;
plotGaussianHotEncodedData(Data, Data_disc, edges, 1);


%% Example Script: Continuous Colored Plot for Data and Data_disc
% Assumptions:
%   - Data is a 3D array with dimensions [time x parameter x trial].
%   - Data_disc is the output from your discretizeAndGaussian_SPparams function.
%       For a discretized parameter p, Data_disc{p} is an array of size 
%       [time x nBins x trial].
%   - For this example we use parameter 1 and trial 1.
%   - current_edges defines the bin edges for parameter 1.
%   - A simple time vector is defined as 1:T.

% %% Settings: Choose Parameter and Trial
paramIndex  = 1;
trialIndex  = 1;
dirIndex=2;

currData=Data;
currData=currData(:,:,:,dirIndex);
currData_disc=Data_disc{paramIndex}(:,:,:,dirIndex);

T = size(currData,1);         % Number of time points
timeVec = (1:T)';        % Time vector (modify if you have a different scale)

% --- Extract the Data trial data (for one parameter and one trial) ---
Data_trial = squeeze(currData(:, paramIndex, trialIndex));  % T x 1

% --- Define discretization edges for parameter 1 (if not already defined) ---
current_edges=edges{paramIndex};
% if ~exist('current_edges','var')
%     % For demonstration: 4 bins spanning the data range (5 edges)
%     current_edges = linspace(min(Data_trial), max(Data_trial), 5);
% end
nBins = length(current_edges) - 1;

% --- Determine bin indices for each time point ---
binIndices = discretize(Data_trial, current_edges);  % T x 1; may contain NaN

% --- Create a colormap with one color per bin ---
colors = lines(nBins);  % nBins x 3 array of distinct RGB colors

% %% Create the Figure with Two Axes
figure;

% %% Top Axes: Plot Data with Continuous, Sloped Segments Colored by Bin
subplot(2,1,1);
hold on;

% For each time point, we plot a segment from (time_i - 0.5) to (time_i + 0.5)
% where the y-values at the boundaries are computed as averages to obtain continuity.
%
% For time point i:
%   - If i > 1, define the left endpoint (x = i - 0.5) as y_left = 0.5*(Data(i-1)+Data(i));
%     otherwise (i==1) use the actual value Data(1).
%   - If i < T, define the right endpoint (x = i + 0.5) as y_right = 0.5*(Data(i)+Data(i+1));
%     otherwise (i==T) use the actual value Data(T).
%
% The segment is drawn in the color corresponding to binIndices(i).

for i = 1:T
    % Compute x-coordinates for the segment endpoints
    if i == 1
        x_left = timeVec(i);      % No previous time point: start at time 1
    else
        x_left = timeVec(i) - 0.5;
    end
    if i == T
        x_right = timeVec(i);     % Last time point: end at time T
    else
        x_right = timeVec(i) + 0.5;
    end
    
    % Compute y-coordinates for the segment endpoints:
    if i == 1
        y_left = Data_trial(i);
    else
        y_left = 0.5 * (Data_trial(i-1) + Data_trial(i));
    end
    if i == T
        y_right = Data_trial(i);
    else
        y_right = 0.5 * (Data_trial(i) + Data_trial(i+1));
    end
    
    % Choose the color based on the bin assigned at time point i.
    currentBin = binIndices(i);
    if isnan(currentBin)
        segColor = [0, 0, 0];    % Use black if no bin was assigned
        segStyle = '--';
    else
        segColor = colors(currentBin, :);
        segStyle = '-';
    end
    
    % Plot the segment
    plot([x_left, x_right], [y_left, y_right], 'Color', segColor, 'LineStyle', segStyle, 'LineWidth', 2);
end

% Optionally, overlay the actual Data values at integer time points for reference
plot(timeVec, Data_trial, 'k:', 'LineWidth', 1);
xlabel('Time');
ylabel('Data Value');
title('Data Time Series (Continuous Segments Colored by Bin)');
hold off;

% %% Bottom Axes: Plot the Corresponding Gaussian Curves from Data_disc
subplot(2,1,2);
hold on;

% --- Extract the discretized (Gaussian) data for this parameter/trial ---
% Here we assume Data_disc{paramIndex} has size [T x nBins x trial]
disc_trial = squeeze(currData_disc(:, :, trialIndex));  % [T x nBins]

% Plot each bin’s Gaussian curve in the same color as used above
for b = 1:nBins
    plot(timeVec, disc_trial(:, b), 'Color', colors(b, :), 'LineWidth', 2);
end

xlabel('Time');
ylabel('Gaussian Amplitude');
title('Data\_disc: Gaussian Curves per Bin');
legend(arrayfun(@(x) sprintf('Bin %d', x), 1:nBins, 'UniformOutput', false), ...
       'Location','best');
hold off;




end
