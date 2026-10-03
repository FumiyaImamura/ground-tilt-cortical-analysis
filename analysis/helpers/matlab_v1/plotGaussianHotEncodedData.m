function plotGaussianHotEncodedData(tpk_numeric, tpk_GaussHot_cellVec, edges_cellVec, paramIndex, options)
arguments
	tpk_numeric			% numeric array
	tpk_GaussHot_cellVec	% cell vector
	edges_cellVec		% cell vector
	paramIndex	% tpk(:,paramIndex,:,...), tpk_disc{paramIndex}, edges{paramIndex} will be used
	options.dirIndex=[]
	options.tiledlayoutHandle=[]
	options.tilenum=[]
	options.Title=[]
end
tpk=tpk_numeric;
tpk_disc=tpk_GaussHot_cellVec;
edges=edges_cellVec;

dirIndex=options.dirIndex;
tiledlayoutHandle=options.tiledlayoutHandle;
tilenum=options.tilenum;
Title=options.Title;
% plotGaussianHotEncodedData - Plot original tpk and Gaussian-hot encoded data with
% continuous colored segments, using vectorized plotting for trial boundaries.
%
% This function produces two subplots:
%
%   1. Top subplot:
%      - The original tpk time series (for a given parameter) is displayed.
%      - All trial/condition dimensions (3rd and higher) are concatenated along the 
%        time axis.
%      - For each time point, a segment is drawn whose endpoints are computed as an
%        average with its neighbors (to yield a continuous slope). Each segment is 
%        colored according to the discretized bin (using edges) for that time point.
%      - Trial boundaries are marked with vertical dashed lines.
%
%   2. Bottom subplot:
%      - The corresponding Gaussian-hot encoded curves from tpk_disc are plotted.
%      - Each bin’s curve is shown in the same color as in the top subplot, and trial
%        boundaries are marked similarly.
%
% Additionally, the x-axes of the two subplots are linked so that you can compare them.
%
% Syntax:
%   plotGaussianHotEncodedData(tpk, tpk_disc, edges, paramIndex, dirIndex)
%
% Inputs:
%   tpk        - Original data, a multi-dimensional array with dimensions:
%                [time x parameter x trial x ...].
%   tpk_disc   - A cell array produced by your discretization function. For each 
%                discretized parameter p, tpk_disc{p} is an array of size 
%                [time x nBins x trial x ...].
%   edges      - A cell array where each cell contains a vector of bin edges for the
%                corresponding parameter.
%   paramIndex - The parameter index to plot (e.g., 1).
%   dirIndex   - (Optional) If your data has a 4th dimension (or higher), use this 
%                index to select a specific slice. If not provided, all data is used.
%
% Example:
%   plotGaussianHotEncodedData(tpk, tpk_disc, edges, 1, 2);

%% see if edges and tpk_disc is cell or not
singular=false;
if ~iscell(tpk_disc)
	if ~iscell(edges)
	singular=true;
	else
		error('edge and tpk_disc data type should match')
	end
else
	if ~iscell(edges)
		error('edge and tpk_disc data type should match')
	end
end

if singular
	paramIndex=1;
	edges={edges};
	tpk_disc={tpk_disc};
end

%% Select the Appropriate Slice if dirIndex is Provided
if ~isempty(dirIndex)
	currtpk = tpk(:,:,:,dirIndex);
	curr_tpk_disc = tpk_disc{paramIndex}(:,:,:,dirIndex);
else
	currtpk = tpk;
	curr_tpk_disc = tpk_disc{paramIndex};
end


%% Concatenate All Trial/Condition Dimensions into the Time Axis
% currtpk is assumed to have dimensions: [T x parameters x trial x ...]
sz = size(currtpk);
T = sz(1);  % number of time points per trial
if numel(sz) > 2
    nTrials = prod(sz(3:end));
else
    nTrials = 1;
end

% Reshape the original data for the specified parameter into [T x nTrials]
data = reshape(currtpk(:, paramIndex, :), [T, nTrials]);  % [T x nTrials]
% Concatenate along time: [T*nTrials x 1]
data_cat = data(:);

% Similarly, reshape the discretized data for the specified parameter.
% tpk_disc dimensions: [T x nBins x trial x ...]
nBins = size(curr_tpk_disc, 2);  % number of bins from tpk_disc output
disc = reshape(curr_tpk_disc, [T, nBins, nTrials]);  % [T x nBins x nTrials]
% Concatenate along time: [T*nTrials x nBins]
disc_cat=BlockCat(disc,3,1);
% % % % % disc_cat = reshape(disc, [T*nTrials, nBins]);

% Retrieve the bin edges for the parameter and update nBins accordingly.
current_edges = edges{paramIndex};
nBins = length(current_edges) - 1;  % number of bins based on edges

% Create a concatenated time vector.
totalTime = T * nTrials;
timeVec = (1:totalTime)';

% Compute the bin index for each time point in the concatenated data.
binIndices = discretize(data_cat, current_edges);

% Create a colormap with one color per bin.
colors = lines(nBins);

%% Vectorized Accumulation of Segments for the Top Subplot
% We accumulate the segment endpoints (x and y) for each bin in cell arrays.
segX = cell(nBins, 1);  % each cell will hold the x coordinates for segments
segY = cell(nBins, 1);  % each cell will hold the corresponding y coordinates

% Reshape binIndices into a [T x nTrials] matrix.
binIndicesMat = reshape(binIndices, [T, nTrials]);

% Loop over each trial (expected to be relatively few compared to total time points)
for trial = 1:nTrials
    % Indices for the current trial in the concatenated vector.
    trial_idx = (trial-1)*T + (1:T)';  % column vector of indices
    
    % Local time and data for this trial.
    localTime = trial_idx;           % global time indices for this trial
    localData = data(:, trial);        % [T x 1]
    
    % Compute x coordinates for segment endpoints:
    % For the first time point, x_left equals localTime(1); for others, subtract 0.5.
    x_left = localTime;
    if T > 1
        x_left(2:end) = localTime(2:end) - 0.5;
    end
    % For the last time point, x_right equals localTime(T); for others, add 0.5.
    x_right = localTime;
    if T > 1
        x_right(1:end-1) = localTime(1:end-1) + 0.5;
    end
    
    % Compute y coordinates for segment endpoints:
    % For the first time point, y_left equals localData(1); for others, average with previous.
    y_left = localData;
    if T > 1
        y_left(2:end) = 0.5 * (localData(1:end-1) + localData(2:end));
    end
    % For the last time point, y_right equals localData(T); for others, average with next.
    y_right = localData;
    if T > 1
        y_right(1:end-1) = 0.5 * (localData(1:end-1) + localData(2:end));
    end
    
    % Get the bin indices for this trial.
    localBins = binIndicesMat(:, trial);
    
    % For each bin, accumulate the segment endpoints.
    for k = 1:nBins
        idx = find(localBins == k);
        if ~isempty(idx)
            % Create segments with NaN between segments for proper plotting.
            segX{k} = [segX{k}; reshape([x_left(idx)'; x_right(idx)'; NaN(1, numel(idx))], [], 1)];
            segY{k} = [segY{k}; reshape([y_left(idx)'; y_right(idx)'; NaN(1, numel(idx))], [], 1)];
        end
    end
end

% %% Create the Figure and Plot the Data
% 
% figure;

%% Top Subplot: Original tpk Data with Continuous Colored Segments
if ~isempty(tiledlayoutHandle)
	tiler=tiledlayout(tiledlayoutHandle,2,1,TileSpacing="tight",Padding="tight");
	tiler.Layout.Tile=tilenum;	
else
	tiler=tiledlayout(2,1,TileSpacing="tight",Padding="tight");
end
tiler.Title.String=Title;
ax1=nexttile(tiler);

hold(ax1, 'on');
% Plot the segments for each bin (one plot call per bin).
for k = 1:nBins
    if ~isempty(segX{k})
        plot(ax1, segX{k}, segY{k}, 'Color', colors(k,:), 'LineWidth', 2);
    end
end
% Overlay the complete concatenated tpk data as a dotted line.
plot(ax1, timeVec, data_cat, 'k:', 'LineWidth', 1);

% ---- Vectorized Plotting of Trial Boundaries for the Top Plot ----
% Build vertical line segments for trial boundaries.
trial_boundaries = (1:(nTrials-1)) * T + 0.5;  % boundaries between trials
y_min_top = min(data_cat);
y_max_top = max(data_cat);
% Create matrices for vertical segments (each column is a segment).
X_bound_top = [trial_boundaries; trial_boundaries; nan(1, numel(trial_boundaries))];
Y_bound_top = [repmat(y_min_top, 1, numel(trial_boundaries));
               repmat(y_max_top, 1, numel(trial_boundaries));
               nan(1, numel(trial_boundaries))];
plot(ax1, X_bound_top(:), Y_bound_top(:), 'k--', 'LineWidth', 1);
% Mark the end of each trial with scatter markers.
trialEndTimes = (T:T:totalTime)';
trialEndValues = data_cat(T:T:totalTime);
scatter(ax1, trialEndTimes, trialEndValues, 40, 'k', 'filled');

xlabel(ax1, 'Time');
ylabel(ax1, 'tpk Value');
title(ax1, 'tpk Time Series (Continuous Segments Colored by Bin)');
hold(ax1, 'off');

%% Bottom Subplot: Gaussian Curves from tpk_disc
ax2=nexttile(tiler);
hold(ax2, 'on');
% Plot each bin's Gaussian curve (data already concatenated) in one call per bin.
for b = 1:nBins
    plot(ax2, timeVec, disc_cat(:, b), 'Color', colors(b,:), 'LineWidth', 2);
end
% ---- Vectorized Plotting of Trial Boundaries for the Gaussian Plot ----
y_min_gauss = min(disc_cat(:));
y_max_gauss = max(disc_cat(:));
X_bound_gauss = [trial_boundaries; trial_boundaries; nan(1, numel(trial_boundaries))];
Y_bound_gauss = [repmat(y_min_gauss, 1, numel(trial_boundaries));
                 repmat(y_max_gauss, 1, numel(trial_boundaries));
                 nan(1, numel(trial_boundaries))];
plot(ax2, X_bound_gauss(:), Y_bound_gauss(:), 'k--', 'LineWidth', 1);

xlabel(ax2, 'Time');
ylabel(ax2, 'Gaussian Amplitude');
title(ax2, 'tpk\_disc: Gaussian Curves per Bin');
legend(ax2, arrayfun(@(x) sprintf('Bin %d', x), 1:nBins, 'UniformOutput', false), 'Location', 'best');
hold(ax2, 'off');

%% Link the x-axes of the two subplots for comparison
linkaxes([ax1, ax2], 'x');

end