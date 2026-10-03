function G=HotEncoder(Data_numeric,encoding_dim,parameter_dim_mustBeSize1,kernel_1d,bincount)
% GaussianHotEncoder_usingDiscretizedData applies a Gaussian smoothing to discretized data.
% It converts the discretized values along the specified parameter dimension 
% into Gaussian curves along the encoding dimension (typically time). The smoothing width is 
% controlled by sigma and the number of discrete bins by bincount.
%
% Input:
%   Data_numeric         - Numeric array with discretized data
%   encoding_dim - Dimension along which the Gaussian is applied
%   parameter_dim_mustBeSize1- Dimension that encoded parameters are going to align in (must be size 1)
%   sigma        - Standard deviation of the Gaussian curve
%   bincount     - Number of bins for the discretization
%
% Output:
%   Data_disc_numeric    - The Gaussian-smoothed data array

arguments
    Data_numeric            % numeric array (discretized). size(Data_numeric,parameter_dim_mustBeSize1) must be 1 since this is the dimension that gaus hot encoding parameters are going to be
    encoding_dim
    parameter_dim_mustBeSize1    
    kernel_1d (1,:)
    bincount
end

Data=Data_numeric;
param_dim=parameter_dim_mustBeSize1;
nBins=bincount;

% Validate
sz = size(Data);
if sz(param_dim)~=1, error('param_dim must be size 1'); end
if encoding_dim==param_dim, error('encoding_dim and param_dim must differ'); end

% Permute to [time, param, rest...]
dims  = 1:ndims(Data);
pvec  = [encoding_dim, param_dim, setdiff(dims,[encoding_dim, param_dim])];
X     = permute(Data, pvec);
rest  = sz(setdiff(dims,[encoding_dim, param_dim]));
T     = sz(encoding_dim);
X     = reshape(X, T, 1, []);       % → [T×1×prod(rest)]

% One‑hot mask [T×B×N]
mask   = X == reshape(1:nBins,1,[],1);	% [T×B×N]

% convolve with kernel
S=FastConv(mask,kernel_1d(:),1,shape="same_as_a");
S=real(S);

% % Build Gaussian kernel (±3σ)
% w = ceil(3*sigma);
% g = exp(-0.5*((-w:w)/sigma).^2)';
% g = g./sum(g);
% 
% % Convolve along time via filter (operates on columns)
% M = reshape(mask, T, []);           % [T×(B·N)]
% S = filter(g, 1, M);
S = reshape(S, T, nBins, []);       % [T×B×N]

% Restore original shape and inverse‑permute
S = reshape(S, [T, nBins, rest]);
G = ipermute(S, pvec);

end