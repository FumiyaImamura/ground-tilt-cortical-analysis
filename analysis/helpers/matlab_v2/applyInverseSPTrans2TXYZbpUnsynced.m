function tXYZbp_out = applyInverseSPTrans2TXYZbpUnsynced(SP_T, tXYZbp, org, opts)
arguments
    SP_T (4,4,:) double
    tXYZbp cell
    org (1,1) struct
    opts.spFirstMasterIndex (1,1) double = 1
    opts.outOfRangePolicy (1,1) string {mustBeMember(opts.outOfRangePolicy, ["nan","nearest","error"])} = "nearest"
end
% APPLYINVERSESPTRANS2TXYZBPUNSYNCED
% Apply inverse platform transform to posture coordinates when SP_T and
% posture samples are not perfectly synchronized.
%
% INPUTS
% ------
% SP_T :
%   4 x 4 x Nsp homogeneous transforms sampled on the 1000 Hz master clock.
%   Assumed to map platform coordinates to lab coordinates:
%
%       p_lab = R * p_platform + t
%
% tXYZbp :
%   cell array; each cell contains a [Nt x 3] position array in lab coordinates.
%
% org.t_fdlc :
%   [Nt x 1] vector. For each posture sample, gives its master-clock index.
%   Values may be non-integer and may be negative.
%
% OPTIONS
% -------
% opts.spFirstMasterIndex :
%   Master-clock index corresponding to SP_T(:,:,1).
%   Default = 0.
%
% opts.outOfRangePolicy :
%   What to do when org.t_fdlc is outside the SP_T range:
%     "nan"     -> output rows become NaN
%     "nearest" -> clamp to nearest SP_T endpoint
%     "error"   -> throw an error
%
% OUTPUT
% ------
% tXYZbp_out :
%   cell array with the same structure as tXYZbp, but positions are expressed
%   in platform coordinates, still sampled at the original posture timestamps.
%
% NOTES
% -----
% - This function uses org.t_fdlc because it directly maps each posture sample
%   to master-clock time.
% - org.fdlc_t is not needed for this operation.
% - Rotation interpolation is done by normalized linear interpolation of
%   sign-corrected unit quaternions, which is appropriate here because SP_T
%   is sampled densely at 1000 Hz.

    if ~isfield(org, "t_fdlc")
        error("org must contain the field org.t_fdlc.");
    end

    tq = org.t_fdlc(:);   % posture-sample times in master-clock index units
    Nt = numel(tq);

    % validateattributes(tq, {"double"}, {"vector","real","finite"}, mfilename, "org.t_fdlc");

	% Validate tXYZbp
	for i = 1:numel(tXYZbp)
		Xi = tXYZbp{i};

		if ~isnumeric(Xi) || size(Xi,2) ~= 3
			error("Each cell of tXYZbp must contain a numeric array with 3 columns.");
		end

		nXi = size(Xi, 1);

		if nXi < Nt
			% org.t_fdlc is longer than this tXYZbp cell: keep Xi as-is
			warning( ...
				"Cell %d of tXYZbp has %d rows, but org.t_fdlc has %d elements. " + ...
				"Using only the first %d elements of org.t_fdlc for this cell; no rows are shed from tXYZbp.", ...
				i, nXi, Nt, nXi);

		elseif nXi > Nt
			% tXYZbp is longer than org.t_fdlc: shed the extra tail rows
			nShed = nXi - Nt;
			warning( ...
				"Cell %d of tXYZbp has %d rows, but org.t_fdlc has %d elements. " + ...
				"Shedding the last %d row(s) from this tXYZbp cell.", ...
				i, nXi, Nt, nShed);

			tXYZbp{i} = Xi(1:Nt, :);
		end
	end

    Nsp = size(SP_T, 3);
    spIdx = opts.spFirstMasterIndex + (0:Nsp-1)';   % master-clock index of each SP_T sample

    % Extract platform translation samples
    tSamples = squeeze(SP_T(1:3,4,:)).';   % [Nsp x 3]

    % Extract platform rotation samples and convert to quaternions
    RSamples = SP_T(1:3,1:3,:);
    qSamples = rotmPagesToQuat(RSamples);   % [Nsp x 4], format [w x y z]

    % Enforce quaternion sign continuity so interpolation is well behaved
    qSamples = enforceQuaternionContinuity(qSamples);

    % Handle out-of-range times
    switch opts.outOfRangePolicy
        case "error"
            if any(tq < spIdx(1) | tq > spIdx(end))
                error("Some org.t_fdlc values lie outside the SP_T sample range.");
            end
            tqUsed = tq;
            valid = true(size(tq));

        case "nearest"
            tqUsed = min(max(tq, spIdx(1)), spIdx(end));
            valid = true(size(tq));

        case "nan"
            valid = (tq >= spIdx(1)) & (tq <= spIdx(end));
            tqUsed = tq(valid);
    end

    % Interpolate translation and rotation at posture timestamps
    tInterp = nan(Nt, 3);
    qInterp = nan(Nt, 4);

    if any(valid)
        tInterp(valid,:) = interp1(spIdx, tSamples, tqUsed, "linear");

        qInterpValid = interp1(spIdx, qSamples, tqUsed, "linear");
        qInterpValid = normalizeRows(qInterpValid);

        qInterp(valid,:) = qInterpValid;
    end

    % Apply inverse transform to each cell
    tXYZbp_out = cell(size(tXYZbp));

    for i = 1:numel(tXYZbp)
        X = tXYZbp{i};          % [Nt x 3], lab coordinates
        Y = nan(size(X), "like", X);

        rowValid = valid & all(~isnan(X), 2);

        if any(rowValid)
            d = X(rowValid,:) - tInterp(rowValid,:);       % subtract interpolated translation
            Y(rowValid,:) = quatRotateInverse(qInterp(rowValid,:), d);  % apply R'
        end

        tXYZbp_out{i} = Y;
    end
end


function q = rotmPagesToQuat(R)
arguments
    R (3,3,:) double
end
% Convert 3x3xN rotation matrices to N x 4 unit quaternions [w x y z].

    N = size(R, 3);

    r11 = reshape(R(1,1,:), N, 1);
    r12 = reshape(R(1,2,:), N, 1);
    r13 = reshape(R(1,3,:), N, 1);
    r21 = reshape(R(2,1,:), N, 1);
    r22 = reshape(R(2,2,:), N, 1);
    r23 = reshape(R(2,3,:), N, 1);
    r31 = reshape(R(3,1,:), N, 1);
    r32 = reshape(R(3,2,:), N, 1);
    r33 = reshape(R(3,3,:), N, 1);

    q = zeros(N, 4, "like", R);
    tr = r11 + r22 + r33;

    idx1 = tr > 0;
    if any(idx1)
        S = 2 * sqrt(1 + tr(idx1));
        q(idx1,1) = 0.25 * S;
        q(idx1,2) = (r32(idx1) - r23(idx1)) ./ S;
        q(idx1,3) = (r13(idx1) - r31(idx1)) ./ S;
        q(idx1,4) = (r21(idx1) - r12(idx1)) ./ S;
    end

    idx2 = ~idx1 & (r11 >= r22) & (r11 >= r33);
    if any(idx2)
        S = 2 * sqrt(1 + r11(idx2) - r22(idx2) - r33(idx2));
        q(idx2,1) = (r32(idx2) - r23(idx2)) ./ S;
        q(idx2,2) = 0.25 * S;
        q(idx2,3) = (r12(idx2) + r21(idx2)) ./ S;
        q(idx2,4) = (r13(idx2) + r31(idx2)) ./ S;
    end

    idx3 = ~idx1 & ~idx2 & (r22 >= r33);
    if any(idx3)
        S = 2 * sqrt(1 - r11(idx3) + r22(idx3) - r33(idx3));
        q(idx3,1) = (r13(idx3) - r31(idx3)) ./ S;
        q(idx3,2) = (r12(idx3) + r21(idx3)) ./ S;
        q(idx3,3) = 0.25 * S;
        q(idx3,4) = (r23(idx3) + r32(idx3)) ./ S;
    end

    idx4 = ~idx1 & ~idx2 & ~idx3;
    if any(idx4)
        S = 2 * sqrt(1 - r11(idx4) - r22(idx4) + r33(idx4));
        q(idx4,1) = (r21(idx4) - r12(idx4)) ./ S;
        q(idx4,2) = (r13(idx4) + r31(idx4)) ./ S;
        q(idx4,3) = (r23(idx4) + r32(idx4)) ./ S;
        q(idx4,4) = 0.25 * S;
    end

    q = normalizeRows(q);
end


function q = enforceQuaternionContinuity(q)
arguments
    q (:,4) double
end
% Flip quaternion signs as needed so adjacent samples stay on the same branch.

    for k = 2:size(q,1)
        if dot(q(k,:), q(k-1,:)) < 0
            q(k,:) = -q(k,:);
        end
    end
end


function q = normalizeRows(q)
arguments
    q (:,:) double
end

    nrm = sqrt(sum(q.^2, 2));
    bad = nrm <= eps(class(q));
    if any(bad)
        error("Encountered a zero-norm quaternion during interpolation.");
    end
    q = q ./ nrm;
end


function vOut = quatRotateInverse(q, v)
arguments
    q (:,4) double
    v (:,3) double
end
% Apply inverse rotation defined by unit quaternion q = [w x y z].
% Equivalent to vOut = R(q)' * v.

    s = q(:,1);
    u = q(:,2:4);

    uv = sum(u .* v, 2);
    uu = sum(u .* u, 2);

    vOut = 2 * uv .* u ...
         + (s.^2 - uu) .* v ...
         - 2 * s .* cross(u, v, 2);
end