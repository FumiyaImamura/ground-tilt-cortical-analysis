function [Trainer, Tester] = splitArrayToGroups(arr, groupnum, dim2split)
% splitArrayGroups splits an array into groups along a specified dimension.
%
%   [Trainer, Tester] = splitArrayGroups(arr, groupnum, dim2split)
%
%   INPUTS:
%       arr       - A numeric or multi-dimensional array.
%       groupnum  - A positive integer specifying the number of groups.
%       dim2split - The dimension along which to perform the split.
%
%   OUTPUTS:
%       Trainer   - A cell array of length groupnum. Each cell contains the 
%                   array with the selected indices (in dimension dim2split) 
%                   removed (i.e., training data).
%       Tester    - A cell array of length groupnum. Each cell contains the 
%                   array with only the selected indices along dim2split (i.e., 
%                   testing data).
%
%   Example:
%       A = rand(10, 5, 8); % Example 3D array.
%       [TrainGroups, TestGroups] = splitArrayGroups(A, 4, 3);
%

    % Validate groupnum
    if ~isscalar(groupnum) || groupnum < 1 || mod(groupnum, 1) ~= 0
        error('groupnum must be a positive integer scalar.');
    end

    % Validate dim2split
    if ~isscalar(dim2split) || mod(dim2split, 1) ~= 0 || dim2split < 1
        error('dim2split must be a positive integer scalar.');
    end

    % Ensure that the dimension to be split is within the array dimensions.
    if dim2split > ndims(arr)
        error('dim2split exceeds the number of dimensions in arr.');
    end

    % Preallocate cell arrays for Trainer and Tester
    Trainer = cell(1, groupnum);
    Tester  = cell(1, groupnum);

    % Determine the size along the dimension we are splitting.
    dimSize = size(arr, dim2split);

    % Loop over each group to split indices
    for i = 1:groupnum
        % Get test indices: every groupnum-th index starting at i.
        testIdx = i:groupnum:dimSize;
        % Training indices are all indices except those chosen for testing.
        trainIdx = setdiff(1:dimSize, testIdx);

        % Create a colon indexing cell array for each dimension.
        ind = repmat({':'}, 1, ndims(arr));

        % Indexing for test data
        ind{dim2split} = testIdx;
        Tester{i} = arr(ind{:});

        % Indexing for training data
        ind{dim2split} = trainIdx;
        Trainer{i} = arr(ind{:});
    end

	if false
		%% Example test script for the splitArrayGroups function

		% Test script for verifying the splitArrayGroups function

		% Create a sample 3D array, for example, a 10 x 5 x 12 array.
		A = rand(10, 5, 12);

		% Set the number of groups to split along the third dimension.
		groupnum = 3;
		dim2split = 3;

		% Call the splitArrayGroups function.
		[Trainer, Tester] = splitArrayToGroups(A, groupnum, dim2split);

		% Display overall settings.
		fprintf('Testing on a %dx%dx%d array, splitting along dimension %d into %d groups.\n\n', ...
			size(A,1), size(A,2), size(A,3), dim2split, groupnum);

		% Loop over each group and print the size of the training and testing arrays.
		for i = 1:groupnum
			trainSize = size(Trainer{i});
			testSize  = size(Tester{i});

			fprintf('Group %d:\n', i);
			fprintf('  Training data size: [%s]\n', num2str(trainSize));
			fprintf('  Testing data size:  [%s]\n\n', num2str(testSize));
		end

		% Additional verification: ensure the indices for training and testing together
		% cover the complete third dimension without overlap.
		totalSlices = zeros(1, groupnum);
		for i = 1:groupnum
			% Number of testing slices for group i.
			totalSlices(i) = size(Tester{i}, dim2split);
		end
		fprintf('Testing: Sum of test slices across groups = %d (should be approximately equal to %d).\n', ...
			sum(totalSlices), size(A,dim2split));



end

end
