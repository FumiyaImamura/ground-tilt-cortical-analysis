function [A_aligned, B_aligned] = ShiftAndGetOverlappingVec(A, B, shift)
% alignment is done to satisfy A(shift+k)<-->B(k)
% when shift is positive, B is pushed forward relative to A, and when it is
% negative, B is pushed back relative to A.

% Check that A and B are the same length.
n = length(A);
if length(B) ~= n
	error('Vectors A and B must have the same length.');
end

% Determine the valid range of indices k for which
% A(shift+k) and B(k) are both defined.
k_start = max(1, 1 - shift);
k_end   = min(n, n - shift);

% Check if there is an overlapping region.
if k_start > k_end
	error('Shift is out of range; no overlapping region exists.');
end

% Extract the aligned segments.
idx = k_start:k_end;
A_aligned = A(shift + idx);
B_aligned = B(idx);
end
