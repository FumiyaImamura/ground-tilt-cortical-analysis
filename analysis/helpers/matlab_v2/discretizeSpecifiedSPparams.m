function [Discretized_t1k_cellVector, DiscretizingFunc_cellVector ]= discretizeSpecifiedSPparams(tpk, edges_cellVector, params2discretize_stringArray4jointDiscretization)
%discretizeSpecifiedSPparams Discretize SP parameters.
%
%   Discretized_tpk = discretizeSpecifiedSPparams(tpk,edges_cellVector,params2discretize_stringArray4jointDiscretization,options)
%
%   tpk is a high‐dimensional numeric array whose 2nd dimension is of size 6,
%   corresponding (in order) to the SP parameters: ["theta","phi","tau","X","Y","Z"].
%
%   edges_cellVector is a 1×6 cell array where each element holds the bin edges for the 
%   corresponding SP parameter.
%
%   params2discretize_stringArray4jointDiscretization is a (repeating) string
%   input. Each element can be a single parameter name (e.g. "theta") or a vector
%   of parameter names (e.g. ["theta","phi"]). In the first case the function
%   discretizes the given parameter using the corresponding edges_cellVector. In the second
%   case the function “joint‐discretizes” the parameters, meaning that it first
%   discretizes each parameter and then combines the resulting indices into a single
%   joint bin (via a sub2ind‐like procedure).
%
%   Note: tpk is high-dimensional so the code uses a general indexing‐cell for the 2nd
%   dimension.
% 
%  To use DiscretizingFunc_cellVector.discretize and
%  DiscretizingFunc_cellVector.undiscretize, just input tpk (numeric array,
%  n x 6 x k x ...)


arguments
    tpk
    edges_cellVector
end
arguments (Repeating)
    params2discretize_stringArray4jointDiscretization (1,:) string {mustBeMember(params2discretize_stringArray4jointDiscretization, ...
                                            ["theta","phi","tau","X","Y","Z"])}
end

edges=edges_cellVector;

% Define the order of SP parameter names.
baseNames = ["theta", "phi", "tau", "X", "Y", "Z"];

% Package the repeating parameter groups into a cell array.
paramGroups = {params2discretize_stringArray4jointDiscretization{:}};
numGroups = numel(paramGroups);

% Preallocate the cell array for discretized results.
Discretized_t1k_cellVector = cell(1, numGroups);
DiscretizingFunc_cellVector=cell(1, numGroups);
for g = 1:numGroups
    % Each group is a row vector of strings (e.g. "theta" or ["theta","phi"])
    group = paramGroups{g};
    
    % Get the corresponding column indices for tpk from baseNames.
    indices = zeros(1, numel(group));
    for k = 1:numel(group)
        idx = find(baseNames == group(k), 1);
        if isempty(idx)
            error("Invalid parameter name: %s", group(k));
        end
        indices(k) = idx;
    end

	% --- Joint discretization for multiple parameters ---
	subs = repmat({':'}, 1, ndims(tpk));
	subs{2}=indices;
	currdata=tpk(subs{:});
	curredges=edges(indices);
	[JointDescretizedData, F] = Discretize_multidimensional(currdata, curredges, 2);
	disc=JointDescretizedData;
    Discretized_t1k_cellVector{g} = disc;

	% Make Discretization function compatible with original data
	newF=F;
	newF.discretize=@(a)F.discretize(a(subs{:}));
	% newF.undiscretize=@(a)F.undiscretize(a(subs{:}));
	newF.undiscretize=@(a)F.undiscretize(a);
	DiscretizingFunc_cellVector{g}=newF;
end


end
