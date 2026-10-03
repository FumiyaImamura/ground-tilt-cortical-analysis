function newx = BlockCat(x, dim2squish, dim2concat)
% BlockCat Concatenates slices of an array along a specified dimension.
%
%   newx = BlockCat(x, dim2squish, dim2concat) extracts slices of the
%   n-dimensional array x along the dimension specified by dim2squish and
%   concatenates these slices along the dimension specified by dim2concat.
%
%   Example:
%       If x is a 4D array and you call:
%           newx = BlockCat(x, 3, 1);
%       then newx is formed by taking slices x(:,:,i,...) for each index i in
%       the 3rd dimension and concatenating them along the 1st dimension.
%
%   Inputs:
%       x            - n-dimensional input array.
%       dim2squish   - the dimension of x to "squish" (i.e., the dimension
%                      whose slices will be extracted).
%       dim2concat   - the dimension along which to concatenate the extracted
%                      slices.
%
%   Output:
%       newx         - the resulting array after concatenation.
%
%   Note:
%       dim2squish must not be equal to dim2concat.


    % Validate inputs: dim2squish should not equal dim2concat.
    if dim2squish == dim2concat
        error('BlockCat:InvalidInput', 'dim2squish should not equal dim2concat');
    end

    % Determine the other dimensions in the array (all dimensions except dim2squish)
    otherDims = setdiff(1:ndims(x), dim2squish);

    % Convert x into a cell array where each cell contains one slice along
    % the dim2squish dimension. The function num2cell splits x along the
    % specified dimensions
    xCells = num2cell(x, otherDims);

    % Concatenate the slices stored in the cell array along the dim2concat
    % dimension. The syntax xCells{:} expands the cell array into a comma-separated
    % list of arrays.
    newx = cat(dim2concat, xCells{:});
    


if false
	%% test this function
	%% Testing BlockCat

	disp('--- Testing BlockCat ---');
	% Create a sample 3D array (2 x 3 x 4)
	x = reshape(1:24, [2, 3, 4]);
	disp('Original array x:');
	disp(x);

	% Use BlockCat to extract slices along dimension 3 and concatenate along dimension 1.
	newx_BlockCat = BlockCat(x, 3, 1);
	disp('Result of BlockCat(x, 3, 1):');
	disp(newx_BlockCat);
	disp(['Size of newx_BlockCat: ' mat2str(size(newx_BlockCat))]);
	%%
	A=rand(2,3,4,5);
	for i=1:ndims(A)
		for j=1:ndims(A)
			if j==i
				continue;
			end
			B=BlockCat(A,j,i);
			C=BlockCatRevert(B,j,i,size(A));
			if ~isequal(A,C)
				error("error at i=%d, j=%d",i,j);
			end
		end
	end
	fprintf("No error!\n")
end

%% Intuitive
% newx=[];
% IND=repmat({':'},1,ndims(x));
% for i=1:size(x,dim2squish)
% 	IND{dim2squish}=i;
% 	newx=cat(dim2concat,newx,x(IND{:}));
% end


%% complex but optimized
% % move dim2concat to the 1st, and dim2squish to 2nd dimension
% % watch out for when dim2squish==1 or dim2concat==1;
% % first, flip 2nd dimension and dim2squish
% pvec1=[1:length(size(x))];
% pvec1(dim2squish)=2;
% pvec1(2)=dim2squish;
% newdim2concat=pvec1(dim2concat);	% make sure to update newdim2concat because pvec1 could have moved dim2concat to dim2squish (only when dim2concat==2)
% 
% % second, flip the 1st dimension and dim2concat
% pvec2=[1:length(size(x))];
% pvec2(newdim2concat)=1;
% pvec2(1)=newdim2concat;
% 
% px1=permute(x,pvec1);
% px2=permute(px1,pvec2);
% px=px2;
% 
% 
% sizepx=size(px);
% newsize=[sizepx(1)*sizepx(2),1,sizepx(3:end)];
% 
% rpx=reshape(px,newsize);
% 
% newx=permute(rpx,pvec2);	% first flip 1st and newdim2concat
% newx=permute(newx,pvec1);	% then flip 2nd (or newdim2concat) and dim2squish

end