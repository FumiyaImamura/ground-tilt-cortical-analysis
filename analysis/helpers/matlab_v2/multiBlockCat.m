function [newx,info]=multiBlockCat(x,dim2squish_vec,dim2concat)
% Ensure that dim2concat and all entries in dim2squish_vec are unique.
combinedDims = [dim2concat, dim2squish_vec(:)'];
if numel(unique(combinedDims)) < numel(combinedDims)
	error('multiBlockCat:InvalidInput', 'dim2squish_vec and dim2concat must be unique.');
end

% dim2squish_vec is a vector of dimensions to squish, in order
dim2squish_vec=dim2squish_vec(:);

newx=x;
sizeEachStep=size(x);
for i=1:length(dim2squish_vec)
	newx=BlockCat(newx,dim2squish_vec(i),dim2concat);
	sizenewx=size(newx);
	currsize=1+zeros(1,size(sizeEachStep,2));
	currsize(1:length(sizenewx))=sizenewx;		% do this so that any trailing dimension sizes are 1
	sizeEachStep=cat(1,sizeEachStep,currsize);
end
info=[];
info.dim2squish_vec=dim2squish_vec;
info.dim2concat=dim2concat;
info.sizeEachStep=sizeEachStep(1:end-1,:);		% save the size of the array each step. Keep in mind you don't need the size of the final array for reconstruction


if false
	%% Testing multiBlockCat

	disp(' ');
	disp('--- Testing multiBlockCat ---');
	% Create a sample 4D array (2 x 3 x 4 x 5)
	x4d = reshape(1:120, [2, 3, 4, 5]);
	disp('Original 4D array x4d:');
	disp(x4d);

	% Call multiBlockCat with dim2squish_vec = [3 2] and concatenate along dimension 1.
	% In the first iteration, slices along dimension 3 will be stacked along dim1;
	% then in the resulting array, slices along dimension 2 (of the current array) are stacked along dim1.
	[newx_multi, info] = multiBlockCat(x4d, [3 2], 1);
	disp('Result of multiBlockCat(x4d, [3 2], 1):');
	disp(newx_multi);

	disp('Intermediate sizes at each step (info.sizeEachStep):');
	disp(info.sizeEachStep);
	disp(['Final size of newx_multi: ' mat2str(size(newx_multi))]);
	%%
	%%
	A=rand(2,3,4,5,6);
	for i=1:ndims(A)
		for j=1:ndims(A)
			if j==i
				continue;
			end
			for k=1:ndims(A)
				if k==i || k==j
					continue;
				end
				for l=1:ndims(A)
					if l==i || l==j || l==k
						continue;
					end
					[B,info]=multiBlockCat(A,[j,k],i);
					C=multiBlockCatRevert(B,info);
					if ~isequal(A,C)
						error("error at i=%d, j=%d",i,j);
					end
				end
			end
		end
	end
	fprintf("No error!\n")
end

end