function GCCA=gcca(X,options)
% gcca  Perform generalized canonical correlation analysis (GCCA)
%
%   gcca(X) computes the MAXVAR formulation of GCCA for a cell array of data
%   blocks X. Each X{i} is assumed to be a matrix of observations (rows)
%   by variables (columns). The function:
%     1) Projects each block onto its null space to build joint structure.
%     2) Solves an eigenproblem to find shared directions Y.
%     3) Computes canonical components A{i} for each block via backslash.
%     4) Projects the data onto these components to obtain variates Z{i}.
%
%   Inputs:
%     X     Cell array of data matrices (n_samples-by-n_features_i).
%
%   Outputs (returned in workspace as GCCA struct):
%     GCCA.Y   Shared basis directions (columns sorted by importance).
%     GCCA.A   Cell array of weight matrices for each group. (canonical components)
%				GCCA.A is made so that GCCA.Z{i}=(X{i}-mean(X{i},1))*GCCA.A{i}
%     GCCA.Z   Cell array of canonical variates for each group. (centered data projected to canonical components)
%	  GCCA.mean:	mean of X
%	  GCCA.std:		std of X

%%%%%%%%%%%%%%%%%
% How to use the data
% GCCA.Z{i}=(X{i}-mean(X{i},1))*GCCA.A{i}
% --> this GCCA.Y can be solved using this GCCA.Z

arguments (Repeating)
	X
end
arguments
	options.NormalizeCanonicalVariates	(1,1) logical =false
	% Make this true to make the output almost exactly the same as canoncorr() (only the sign of the vector may be different). 
	% Keep this false if you intend on keeping the alignment
	% between Z and Y, since Y is derived from unnormalized Z, and
	% normalizing Z columnwise across all the groups leads to recalculation of Y
	% which has worse correlation.
end

p_i = cellfun(@(Xi) size(Xi,2), X);		% number of variables per block
K_max = min( size(X{1},1), min(p_i) );	% the absolute upper‐bound on components

% Keep raw data for later projection
rawX = X;

% Compute means and standard deviations, then z-score each block
mu    = cellfun(@(a) mean(a,1), rawX, UniformOutput=false);
sigma = cellfun(@(a) std(a,0,1), rawX, UniformOutput=false);
X    = cellfun(@(a,m,s) (a - m) ./ s, rawX, mu, sigma, UniformOutput=false);

% normalize X data
X=cellfun(@(a)zscore(a,0,1),X,UniformOutput=false);

% obtain group configuration (anchor)
P=cellfun(@(X)eye(size(X,1))-orth(X)*orth(X)',X,UniformOutput=false);
P=sum(cat(3,P{:}),3);
P=(P+P')/2;		% enforce exact symmetry

% solve P'y=lambda*y and obtain y with the smallest eigenvalues
[V,D]=eigs(P',K_max,'smallestabs');
d=diag(D);
[~,ind]=sort(d,"ascend");
Y=V(:,ind);
Y=Y(:,1:K_max);	% limit the # of components to only the valid ones

Az=cellfun(@(Xi)Xi\Y,X,UniformOutput=false); % canonical components (Axes for each group)
A = cellfun(@(Ai, s) Ai ./ s', Az, sigma, UniformOutput=false);	% Transform components back to original data scale--> rawX{i}*A{i} yields canonical variants
Z = cellfun(@(Xi, Ai) (Xi-mean(Xi,1)) * Ai, rawX, A, UniformOutput=false);	% Project raw data onto scaled components

if options.NormalizeCanonicalVariates
	% Scale axes so that Z have unit variance
	for i = 1:numel(Z)
		Zi   = Z{i};
		sZi  = std(Zi, 0, 1);
		sZi(sZi == 0) = 1;
		A{i} = A{i} ./ sZi;
		Z{i} = Zi   ./ sZi;
	end
end

% obtain canonical correlations (R{i} == columnwise_correlation( Y, Z{i} ))
R = cellfun(@(Zi) diag(corr(Y, Zi))', Z, 'UniformOutput', false);

GCCA=[];
GCCA.Y=Y;
GCCA.A=A;
GCCA.Z=Z;
GCCA.R=R;
GCCA.mean=mu;
GCCA.std=sigma;

if false
	%% unit test
	a=rand(10,5).*[1:5]+[1:5];
	b=rand(10,5).*[1:5]+[1:5];
	% a=zscore(a,0,1);
	% b=zscore(b,0,1);
	[A,B,r,U,V]=canoncorr(a,b);
	aa=zscore(a,0,1);
	bb=zscore(b,0,1);
	[C,D]=canoncorr(aa,bb);
	GCCA=gcca(a,b);
	altA=GCCA.A{1};
	ratio=A(:,1)./altA(:,1)	% the axes obtained using gcca method is supposed to be parallel to those of standard cca
	ratio=A(:,2)./altA(:,2)
	ratio=A(:,3)./altA(:,3)

	% reconstruct data
	rec=a*A;
	altrec=(a-mean(a,1))*GCCA.A{1};
	isequal(GCCA.Z{1},altrec)

	% U(:,1)./GCCA.Z{1}(:,1)
end


end
 