function DLC_out = applyInverseSPTrans2DLC(DLC, T)
% APPLYINVERSETRANSFORMATION Inverse-transform DLC 3D coordinate data.
%
%   DLC_out = applyInverseTransformation(DLC, T) applies an inverse
%   transformation to the 3D coordinates stored in the structure DLC.
%
%   The input structure DLC is expected to contain three fields:
%       - tXYZbp   : A cell array of position data. Each cell contains an
%                    array of size [N, 3, dim3, dim4, ...] where the second
%                    dimension holds the 3D (x,y,z) coordinates.
%       - tDXYZbp  : A cell array of velocity data, with the same array shape.
%       - tD2XYZbp : A cell array of acceleration data, with the same array shape.
%
%   The transformation T is a homogeneous transformation matrix with size
%       [4, 4, N, dim3, dim4, ...]
%   where the first two dimensions contain the 4-by-4 transform, and the extra
%   (page) dimensions (N, dim3, dim4, ...) correspond to those in each cell
%   array element of DLC.
%
%   For the position data in tXYZbp, the function applies the full inverse:
%
%       p_out = R' * (p_in - t)
%
%   where R is the rotation matrix (the top-left 3-by-3 part of T) and t is
%   the translation vector (the first 3 elements of the 4th column of T).
%
%   For velocity and acceleration data (tDXYZbp and tD2XYZbp), only the rotation
%   is applied:
%
%       v_out = R' * v_in
%
%   The function returns DLC_out, a structure with the same fields as DLC, with
%   the inverse-transformed data.
%
%   Example:
%       DLC_transformed = applyInverseTransformation(DLC, T);
%
%   See also: cellfun, pagemtimes, permute, ipermute.

    %% Dynamically extract the rotation and translation parts from T
    nDimsT = ndims(T);
    idx_all = repmat({':'}, 1, nDimsT);
    
    % Extract the rotation part R: the top-left 3-by-3 block.
    idx_R = idx_all;
    idx_R{1} = 1:3;
    idx_R{2} = 1:3;
    R = T(idx_R{:});
    
    % For an orthogonal rotation matrix, the inverse is its transpose.
    % Permute dimensions so that the first two (3x3) are swapped.
    Rinv = permute(R, [2, 1, 3:nDimsT]);
    
    % Extract the translation vector: first 3 elements of the 4th column.
    idx_t = idx_all;
    idx_t{1} = 1:3;
    idx_t{2} = 4;
    tvec = T(idx_t{:});
    
	%% Apply the inverse transformation to each field in DLC
	% For positions: subtract translation then apply rotation.
	% For velocities/accelerations: apply rotation only.
	DLC_out = DLC;
	DLC_out.tXYZbp = cellfun(@(X) transformData(X, Rinv, tvec), ...
		DLC.tXYZbp, 'UniformOutput', false);
	DLC_out.tDXYZbp = cellfun(@(X) transformData(X, Rinv, []), ...
		DLC.tDXYZbp, 'UniformOutput', false);
	DLC_out.tD2XYZbp = cellfun(@(X) transformData(X, Rinv, []), ...
		DLC.tD2XYZbp, 'UniformOutput', false);
end

%% ------------------------------------------------------------------------
function Y = transformData(X, Rinv, tvec)
% TRANSFORMDATA Apply the inverse transformation to a single data array.
%
%   Y = transformData(X, Rinv, tvec) applies an inverse transformation to the
%   array X. X is assumed to have size [N, 3, dim3, dim4, ...] where the 3
%   in the second dimension represents the 3D coordinates.
%
%   Rinv is the inverse rotation (i.e. the transpose of the rotation matrix)
%   with size [3, 3, N, dim3, ...]. tvec is the translation vector with size
%   [3, 1, N, dim3, ...]. When tvec is nonempty, the function subtracts tvec
%   from X before applying the rotation. When tvec is empty (as in the case
%   for velocity or acceleration data), only the rotation is applied.
%
%   The output Y has the same size as X.

    dims = size(X);
    nDimsX = ndims(X);
    % Compute the number of pages: product of the 1st dimension and dims 3:end.
    nPages = prod(dims([1, 3:nDimsX]));
    
    % Permute and reshape X so that the coordinate (3) dimension comes first.
    % Xp is then 3-by-(#pages).
    Xp = reshape(permute(X, [2, 1, 3:nDimsX]), 3, nPages);
    
    if ~isempty(tvec)
        % Reshape tvec to match the pages.
        tPages = reshape(tvec, 3, nPages);
        diffVec = Xp - tPages;
    else
        diffVec = Xp;
    end
    
    % Reshape Rinv into pages, each of size 3-by-3.
    Rinv_pages = reshape(Rinv, 3, 3, nPages);
    
    % Use pagemtimes to multiply each page:
    % Multiply Rinv_pages (3x3) with diffVec (3x1) for each page.
    Yp = pagemtimes(Rinv_pages, reshape(diffVec, 3, 1, nPages));  % result: 3-by-1-by-nPages
    
    % Reshape Yp back into the original dimensions:
    % First, put the 3 back as the first dimension, then recover dims(1) and the pages.
    Yp = reshape(Yp, [3, dims(1), dims(3:end)]);
    
    % Inverse permute to restore the original order ([N, 3, dim3, ...]).
    Y = ipermute(Yp, [2, 1, 3:nDimsX]);
end
