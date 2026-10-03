function B = ridge_optimized(yraw, Xraw, k, option)
% RIDGEOPTIMIZED Ridge regression with multiple ridge parameters.
%
%   B = ridgeoption(yraw, Xraw, k, option) computes ridge regression
%   coefficients for one or more ridge penalties, k. When k is a scalar,
%   the method used is determined by the input "option" (1, 2, 3, or 4).
%
%   When k is a vector (of length L), the function computes the ridge
%   regression coefficients for each k simultaneously using an eigen-
%   decomposition. The output B is a (p+1) x s x L array, where:
%       - B(1,:,i) is the intercept for k(i),
%       - B(2:end,:,i) are the slope coefficients for k(i).
%
%   If you want to reconstruct a data
%   y_recon=B(1,:,:)+pagemtimes(Xraw,B(2:end,:,:));
%   y_recon(:,:,i) contains reconstructed result with k(i)
% 
%
%   Inputs:
%       yraw  - (n x s) matrix of responses (each column a different response)
%       Xraw  - (n x p) matrix of predictors
%       k     - ridge penalty (scalar or vector)
%       option- (only used if k is scalar) selects which solver to use:
%                   1: MATLAB backslash
%                   2: linsolve with SYM flag
%                   3: Cholesky factorization with two triangular solves
%                   4: linsolve on Cholesky factors with UT/LT flags
%
%   Output:
%       B     - ((p+1) x s x length(k)) array of regression coefficients.
%               For each response:  B(1,:,:) is the intercept and
%               B(2:end,:,:) are the slopes.
%


% --- Data Preprocessing: Center and scale X and center y ---
Xmean = mean(Xraw, 1);
Xstd  = std(Xraw, 0, 1);
ymean = mean(yraw, 1);

% find 0 std parameters
zeroinds=find(Xstd==0);
if ~isempty(zeroinds)
	warning(sprintf("Constant predictor: %s, setting std to 0 and weight 0",strjoin(string(zeroinds),", ")));
	Xstd(zeroinds)=1;
end

% Center y
ycentered = yraw - ymean;
% Scale X (same as: X = zscore(Xraw,0,1))
X = (Xraw - Xmean) ./ Xstd;

% Get dimensions
[n, p] = size(X);
[~, s] = size(yraw);

% --- Precompute cross-products ---
XtX = X' * X;    % (p x p)
Xty = X' * ycentered;  % (p x s)

% --- Case 1: Single k value ---
if numel(k) == 1
    XtX_reg = XtX + k * eye(p);
    switch option
        case 1
            B0 = XtX_reg \ Xty;
        case 2
            opts.SYM = true;
            B0 = linsolve(XtX_reg, Xty, opts);
        case 3
            R = chol(XtX_reg); % Cholesky factor (upper triangular)
            B0 = R \ (R' \ Xty);
        case 4
            R = chol(XtX_reg);
            opts1.UT = true; opts2.LT = true;
            B0 = linsolve(R, linsolve(R', Xty, opts2), opts1);
        otherwise
            error('Option must be 1, 2, 3, or 4.');
    end
    
    % Rescale the slopes back to the original X scale:
    B1 = B0 ./ Xstd';
    % Compute the intercept:
    Bconst = ymean - Xmean * B1;
    % Combine into final coefficient matrix:
    B = [Bconst; B1];
    % Reshape to a 3D array with singleton third dimension.
    B = reshape(B, p+1, s, 1);
    
% --- Case 2: Multiple k values ---
else
	L = numel(k);
	% For multiple k's, use the eigen-decomposition of XtX:
	[Q, Lambda] = eig(XtX);
	lambda = diag(Lambda);  % (p x 1)

	% Compute Q'*Xty once (p x s)
	QXty = Q' * Xty;

	% For each k, we need the diagonal of (lambda + k). Form a (p x L) matrix:
	invDiag = 1 ./ (lambda + reshape(k, 1, L));  % (p x L)

	% We now want to compute for each k:
	%   B0(:,:,i) = Q * ( diag(1./(lambda + k(i))) * (QXty) )
	%
	% To vectorize, we first replicate QXty to a 3D array:
	QXty_rep = repmat(QXty, [1, 1, L]);   % (p x s x L)
	% Reshape invDiag to (p x 1 x L) so it multiplies each row of QXty:
	invDiag = reshape(invDiag, [p, 1, L]); % (p x 1 x L)
	% Multiply element‐wise:
	middle = QXty_rep .* invDiag;         % (p x s x L)

	% Now compute Q * middle for each k. If your MATLAB version supports
	% pagemtimes (R2021b and later), use it; otherwise, loop over pages.

	B0 = pagemtimes(Q, middle);  % (p x s x L)


	% Rescale slopes to the original scale:
	B1 = B0 ./ reshape(Xstd, [p, 1, 1]);  % (p x s x L)

	% Compute intercept for each k:
	% Note: Xmean is 1 x p and B1 is p x s x L.

	Bconst = ymean - pagemtimes(reshape(Xmean, [1, p]), B1);  % (1 x s x L)


	% Stack intercept and slopes so that B is ((p+1) x s x L)
	B = cat(1, Bconst, B1);
end
end