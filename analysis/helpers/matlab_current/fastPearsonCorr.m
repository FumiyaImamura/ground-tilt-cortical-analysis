function corrCoef = fastPearsonCorr(A, B, dim)
% fastPearsonCorr Compute Pearson correlation along a specified dimension.
%
%   corrCoef = fastPearsonCorr(A, B, dim) computes the Pearson correlation
%   coefficients between arrays A and B along the dimension specified by 'dim'.
%
%   The function supports high-dimensional arrays and broadcasting. Two arrays
%   A and B are considered broadcast-compatible if, for each dimension, their sizes
%   are equal or one of them is 1. However, if the broadcasting occurs along the
%   correlation dimension (i.e. one input is singleton along 'dim' while the other
%   has a larger size), an error is thrown.
%
%   Inputs:
%       A   - Numeric array.
%       B   - Numeric array that is either the same size as A or broadcast-compatible.
%       dim - A positive integer specifying the dimension along which to compute the
%             correlation.
%
%   Output:
%       corrCoef - An array of Pearson correlation coefficients computed along 'dim'.
%                  Its size is the broadcasted size of A and B with the size along
%                  dimension 'dim' equal to 1.
%
%   Example:
%       % A and B are 3-D arrays of size 10 x 5 x 3. Correlation is computed along
%       % the first dimension.
%       A = randn(10, 5, 3);
%       B = randn(10, 5, 3);
%       corrCoef = fastPearsonCorr(A, B, 1);
%
%   See also mean, sum, sqrt.

    % Validate number of inputs.
    if nargin < 3
        error('fastPearsonCorr requires three inputs: A, B, and dim.');
    end

    % Validate that A and B are numeric.
    if ~isnumeric(A) || ~isnumeric(B)
        error('A and B must be numeric arrays.');
    end

    % Validate 'dim'
    if ~isscalar(dim) || dim < 1 || floor(dim) ~= dim
        error('dim must be a positive integer scalar.');
    end

    % Determine the maximum number of dimensions to consider.
    nDim = max([ndims(A), ndims(B), dim]);

    % Get sizes of A and B and pad with ones to length nDim.
    szA = size(A);
    szB = size(B);
    if numel(szA) < nDim
        szA(end+1:nDim) = 1;
    end
    if numel(szB) < nDim
        szB(end+1:nDim) = 1;
    end

    % Check broadcasting compatibility for all dimensions.
    for i = 1:nDim
        if ~(szA(i) == szB(i) || szA(i) == 1 || szB(i) == 1)
            error('A and B are not broadcast-compatible.');
        end
    end

    % Check if broadcasting occurs along the correlation dimension.
    if (szA(dim) == 1 && szB(dim) > 1) || (szB(dim) == 1 && szA(dim) > 1)
        error(['Broadcasting along the correlation dimension (dim = %d) ', ...
               'is not allowed.'], dim);
    end

    % Compute means along the specified dimension.
    meanA = mean(A, dim);
    meanB = mean(B, dim);

    % Center the arrays by subtracting the mean (implicit expansion applies).
    A0 = A - meanA;
    B0 = B - meanB;

    % Compute the numerator: sum along 'dim' of the elementwise product.
    numerator = sum(A0 .* B0, dim);

    % Compute the denominator: product of the Euclidean norms along 'dim'.
    denom = sqrt( sum(A0.^2, dim) .* sum(B0.^2, dim) );

    % Compute Pearson correlation coefficients.
    % Where the denominator is zero (e.g., constant data), return NaN.
    corrCoef = numerator ./ denom;
    corrCoef(denom == 0) = NaN;
end
