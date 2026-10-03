function plotSimilarActivity(numericData, options)
arguments
	numericData
	options.offsetParam=[]
	options.offset_prctile=[]
	options.showTrialBoundary=true
	options.labels=[]
	options.labelsize=[]
	options.EnableColor (1,1) logical =true
	options.normalizeData (1,1) logical =false
end
% plotCellActivity - Plot concatenated cell activity across trials.
%
%   plotCellActivity(neuralData, offsetParam) takes a multi-dimensional array
%   'numericData' where the 1st dimension represents time, the 2nd represents
%   cells or parameters, and the 3rd (and higher) dimensions represent trials. It concatenates
%   all trials along the time axis, plots each cell's activity on the same axes
%   with a vertical offset, and marks the beginning of each trial with a vertical
%   dashed line.
%
%   offsetParam is an optional scalar parameter that determines the vertical offset
%   between each cell's activity trace. Smaller values lead to less separation (and
%   potential overlap) while larger values increase separation.
%
%   Example:
%       % Create a sample data set with 100 time points, 10 cells, 5 trials:
%       neuralData = randn(100, 10, 5);
%       plotCellActivity(neuralData, 20);  % Here, an offset of 20 is used.
%

offsetParam=options.offsetParam;
showTrialBoundary=options.showTrialBoundary;
cellLabels=options.labels;
labelsize=options.labelsize;

% Get size information
dims = size(numericData);
T = dims(1);              % number of time points per trial
numCells = dims(2);       % number of cells
numTrials = prod(dims(3:end)); % total number of trials

% Reshape neuralData to a 3D array: time x cells x trials
numericData = reshape(numericData, T, numCells, numTrials);

% Preallocate the concatenated data array: (T * numTrials) x numCells
concatenatedData = zeros(T * numTrials, numCells);
trialStarts = zeros(numTrials, 1);  % to record start index of each trial

% Concatenate trials along the time dimension and record trial boundaries
for trial = 1:numTrials
    idx_start = (trial-1) * T + 1;
    idx_end = trial * T;
    concatenatedData(idx_start:idx_end, :) = numericData(:, :, trial);
    trialStarts(trial) = idx_start; % save trial start index
end

if options.normalizeData
	concatenatedData=(concatenatedData-mean(concatenatedData,1))./std(concatenatedData,0,1);
end

% % Create the concatenated time axis
% timeAxis = (1:(T * numTrials))';
% Sampling frequency is 10 Hz, convert sample indices to seconds.
fs = 10; % Hz
timeAxis = ((0:(T * numTrials - 1))') / fs;

% Determine offset factor
if isempty(offsetParam)
    % Default offset factor set relative to the data range
    offsetFactor = max(abs(concatenatedData(:)));
else
    offsetFactor = offsetParam;
end
if ~isempty(options.offset_prctile)
	offsetFactor=prctile(abs(concatenatedData(:)),options.offset_prctile);
end

if offsetFactor == 0
    offsetFactor = 1; % fallback if data is all zeros
end
offsets = (0:(numCells-1)) * offsetFactor;
offsets=offsets(end:-1:1);

% Plot the cell activities with vertical offsets

for cellIdx = 1:numCells
	if ~options.EnableColor
		plot(timeAxis, concatenatedData(:, cellIdx) + offsets(cellIdx), 'LineWidth', 1.2,'Color',[1,1,1]*0.3);
	else
		plot(timeAxis, concatenatedData(:, cellIdx) + offsets(cellIdx), 'LineWidth', 1.2);
	end
	hold on
end

% Determine y-axis limits (taking offsets into account) for the boundary lines
allY = concatenatedData + repmat(offsets, size(concatenatedData,1), 1);
yLimits = [min(allY(:)), max(allY(:))];

if showTrialBoundary
	% Plot vertical dashed lines to mark the beginning of each trial (except the first)
	% Convert trial start sample indices to seconds.
	trialBoundaryX = ((trialStarts(2:end) - 1) / fs);
	numBoundaries = length(trialBoundaryX);

	% Prepare matrices with NaNs inserted between segments to avoid unwanted connections
	Xlines = nan(3, numBoundaries);
	Ylines = nan(3, numBoundaries);

	Xlines(1,:) = trialBoundaryX;
	Xlines(2,:) = trialBoundaryX;
	Ylines(1,:) = yLimits(1);
	Ylines(2,:) = yLimits(2);
	% Third row remains NaN to break the connection between segments

	plot(Xlines(:), Ylines(:), 'k--', 'LineWidth', 1);
end

% Process cell labels;
if ~isempty(cellLabels)
	if length(cellLabels) ~= numCells
		error('Number of cell labels must match number of cells.');
	end

	% Determine a suitable x position for the labels (slightly left of the first time point)
	xLabelPos = timeAxis(1) - 0.05 * (timeAxis(end) - timeAxis(1));

	textHandles = cell(numCells, 1);
	% Add text labels for each cell; align vertically with the mean of the cell's data.
	for cellIdx = 1:numCells
		yLabelPos = offsets(cellIdx) + mean(concatenatedData(:, cellIdx));
		if isempty(labelsize)
			textHandles{cellIdx} =text(xLabelPos, yLabelPos, cellLabels(cellIdx), 'HorizontalAlignment', 'right', ...
				'VerticalAlignment', 'middle', 'FontSize', 10);
		else
			textHandles{cellIdx} =text(xLabelPos, yLabelPos, cellLabels(cellIdx), 'HorizontalAlignment', 'right', ...
			'VerticalAlignment', 'middle', 'FontSize', labelsize);
		end
	end
end

% Force a draw to update text extents.
drawnow;

% Adjust x-axis limits to ensure labels are fully visible.
minLeft = inf;
for cellIdx = 1:numCells
    extent = get(textHandles{cellIdx}, 'Extent');  % extent = [x, y, width, height]
    minLeft = min(minLeft, extent(1));
end

% Add some padding based on the total time range.
xPadding = 0.05 * (timeAxis(end) - timeAxis(1));
xlim([minLeft - xPadding, timeAxis(end)]);

% Set y-axis limits so that the data exactly spans the entire axes.
ylim(yLimits+[-offsetFactor,offsetFactor]);

xlabel('Time (s)');
% ylabel('Cell Activity');
% title('Cell Activity Across Trials (concatenated)');
hold off;
set(gca,"YTick",[])

end
