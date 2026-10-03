function newx=InsertSingleDim(x,n)
% This function inserts a singleton dimension to x in a dimension specified
% by n, and shifts the remaining dimensions backwards
dim=n;

sizex=1+zeros(1,max(length(size(x)),dim));	% if dim is greater than non-1 x dim, this max func is necessary
sizex(1:length(size(x)))=size(x);

permvec=[1:length(sizex)];
pre_n=permvec(1:n-1);
post_n=permvec(n:end);
permvec=[pre_n,length(sizex)+1,post_n];
newx=permute(x,permvec);

end