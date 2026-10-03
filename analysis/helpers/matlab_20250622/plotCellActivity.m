function plotCellActivity(neuralData, offsetParam,showTrialBoundary)
arguments
	neuralData
	offsetParam=[]
	showTrialBoundary=true
end
% plotCellActivity - Plot concatenated cell activity across trials.
%
%   plotCellActivity(neuralData, offsetParam) takes a multi-dimensional array
%   'neuralData' where the 1st dimension represents time, the 2nd represents
%   cells, and the 3rd (and higher) dimensions represent trials. It concatenates
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

% Get size information
dims = size(neuralData);
T = dims(1);              % number of time points per trial
numCells = dims(2);       % number of cells
numTrials = prod(dims(3:end)); % total number of trials

% Reshape neuralData to a 3D array: time x cells x trials
neuralData = reshape(neuralData, T, numCells, numTrials);

% Preallocate the concatenated data array: (T * numTrials) x numCells
concatenatedData = zeros(T * numTrials, numCells);
trialStarts = zeros(numTrials, 1);  % to record start index of each trial

% Concatenate trials along the time dimension and record trial boundaries
for trial = 1:numTrials
    idx_start = (trial-1) * T + 1;
    idx_end = trial * T;
    concatenatedData(idx_start:idx_end, :) = neuralData(:, :, trial);
    trialStarts(trial) = idx_start; % save trial start index
end


% % Create the concatenated time axis
% timeAxis = (1:(T * numTrials))';
% Sampling frequency is 10 Hz, convert sample indices to seconds.
fs = 10; % Hz
timeAxis = ((0:(T * numTrials - 1))') / fs;

% Determine offset factor
if isempty(offsetParam)
    % Default offset factor set relative to the data range
    offsetFactor = max(abs(concatenatedData(:))) * 1.5;
else
    offsetFactor = offsetParam;
end
if offsetFactor == 0
    offsetFactor = 1; % fallback if data is all zeros
end
offsets = (0:(numCells-1)) * offsetFactor;

% Plot the cell activities with vertical offsets

for cellIdx = 1:numCells
    plot(timeAxis, concatenatedData(:, cellIdx) + offsets(cellIdx), 'LineWidth', 1.2,'Color',[1,1,1]*0.3);
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

xlabel('Time (s)');
ylabel('Cell Activity');
title('Cell Activity Across Trials (concatenated)');
hold off;
set(gca,"YTick",[])

end
