function orgx=BlockCatRevert(x,dim2squish,dim2concat,OriginalSize)
% Undo concatenation by BlockCat. set the same dim2squish and dim2concat
% you used for BlockCat. The squished dimension should have the size of 1. 
% OriginalSize is the size of the original data before applying BlockCat.


% move dim2concat to the 1st, and dim2squish to 2nd dimension

% permvec=[1:length(size(x))];
% permvec(1)=dim2concat;
% permvec(dim2concat)=1;
% permvec(2)=dim2squish;
% permvec(dim2squish)=2;


% watch out for when dim2squish==1 or dim2concat==1;
% first, flip 2nd dimension and dim2squish
numdims=max([numel(OriginalSize),dim2squish,dim2concat]);
OriginalSize(end+1:numdims)=1;

pvec1=[1:numdims];
pvec1(dim2squish)=2;
pvec1(2)=dim2squish;
newdim2concat=pvec1(dim2concat);	% make sure to update newdim2concat because pvec1 could have moved dim2concat to dim2squish (only when dim2concat==2)

% second, flip the 1st dimension and dim2concat
pvec2=[1:numdims];
pvec2(newdim2concat)=1;
pvec2(1)=newdim2concat;

px1=permute(x,pvec1);
px2=permute(px1,pvec2);
px=px2;

% px=permute(x,permvec);

sizepx=size(px);
newsize=[OriginalSize(dim2concat),OriginalSize(dim2squish),sizepx(3:end)];

rpx=reshape(px,newsize);

% orgx=permute(rpx,permvec);
newx=permute(rpx,pvec2);	% first flip 1st and newdim2concat
orgx=permute(newx,pvec1);	% then flip 2nd (or newdim2concat) and dim2squish

if false
	%% test code
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
end
end
