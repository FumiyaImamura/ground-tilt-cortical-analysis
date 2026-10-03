function newx=DeleteSingleDim(x,dim)
% This function deletes a singleton dimension to x in a dimension specified
% by dim, and shifts the remaining dimensions backwards
n=dim;

if size(x,n)>1
	disp('Error: The size of dim is greater than 1 !!!');
end

pre_n=[1:n-1];
post_n=[n+1:ndims(x)];
permvec=[pre_n,post_n,n];
newx=permute(x,permvec);

end