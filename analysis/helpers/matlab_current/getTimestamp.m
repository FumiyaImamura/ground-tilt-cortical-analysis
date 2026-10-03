function ts = getTimestamp(includeSeconds)
% genTimestamp  Generate a date-time string for file names.
%
%   ts = genTimestamp() returns a string 'yyyy-MM-dd_HH-mm' using the
%   current date and time (hours and minutes), e.g. '2025-07-15_14-32'.
%
%   ts = genTimestamp(true) returns 'yyyy-MM-dd_HH-mm-ss', e.g.
%   '2025-07-15_14-32-37'.  Seconds are omitted if includeSeconds=false.
%
%   Example:
%     >> genTimestamp(true)
%     ans =
%       '2025-07-15_14-32-37'

    if nargin < 1
        includeSeconds = false;
    end

    if includeSeconds
        fmt = 'yyyy-MM-dd_HH-mm-ss';
    else
        fmt = 'yyyy-MM-dd_HH-mm';
    end

    ts = string(datetime('now', 'Format', fmt));
end
