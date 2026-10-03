function vec = force_monotonic_increase(vec_in, options)
arguments
    vec_in (:,1) double
    options.max (1,1) double = Inf
    options.min (1,1) double = -Inf
	options.allow_flat=false
end

% Preserve original orientation
wasRow = isrow(vec_in);
vec = vec_in(:);

n = numel(vec);
if n < 2
    if wasRow
        vec = vec.';
    end
    return;
end

% Compute differences and monotonic mask
d = diff(vec);
if options.allow_flat
	mask = d >= 0;
else
	mask = d > 0;
end


% Identify contiguous monotonic regions
runs = struct('start',{},'length',{});
i = 1;
while i <= numel(mask)
    if mask(i)
        s = i; l = 1;
        i = i + 1;
        while i <= numel(mask) && mask(i)
            l = l + 1;
            i = i + 1;
        end
        runs(end+1).start = s; %#ok<*AGROW>
        runs(end).length = l;
    else
        i = i + 1;
    end
end

if ~isempty(runs)
    % Use the longest monotonic segment
    [~, idx] = max([runs.length]);
    S = runs(idx).start;
    L = runs(idx).length;
    E = S + L - 1;

    % Adjust head if needed
    if S > 1
        slope = d(S);
        for k = S:-1:1
            vec(k) = vec(k+1) - slope;
        end
    end

    % Adjust tail if needed
    if E < n-1
        slopeTail = d(E);
        for k = (E+2):n
            vec(k) = vec(k-1) + slopeTail;
        end
    end
end

% Cap values to bounds
vec(vec > options.max) = options.max;
vec(vec < options.min) = options.min;

% Restore orientation
if wasRow
    vec = vec.';
end

end
