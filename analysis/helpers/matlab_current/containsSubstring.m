% CONTAINSSUBSTRING Determines if any of the substrings in str2search are found in each element of strArray.
%   tf = containsSubstring(strArray, str2search) returns a logical array of the same size as strArray,
%   where each element is true if any substring from str2search is contained in the corresponding element of strArray.
%   The optional name-value pair 'CapSensitive' controls case sensitivity (default is true for case sensitive search).

function tf = containsSubstring(strArray, str2search, options)
arguments
    strArray string
    str2search string
    options.CapSensitive = true
end

if options.CapSensitive
    ignore = false;
else
    ignore = true;
end

% Initialize the output logical array.
tf = false(size(strArray));

% Loop over each search substring and combine results.
for k = 1:numel(str2search)
    tf = tf | contains(strArray, str2search(k), 'IgnoreCase', ignore);
end

if false
	%% test code
	% Example to test containsSubstring function

	strArray = ["apple pie", "banana bread", "cherry tart", "grape juice"];
	str2search = ["apple", "banana", "orange"];

	% Case sensitive search (default)
	tf1 = containsSubstring(strArray, str2search);
	disp(tf1)
	% Expected output: [true  true  false  false]

	% Case insensitive search
	tf2 = containsSubstring(strArray, ["APPLE", "BANANA"], CapSensitive=false);
	disp(tf2)
	% Expected output: [true  true  false  false]


end

end
