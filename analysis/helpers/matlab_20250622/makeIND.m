function IND=makeIND(numericdata,dim,ind4EachDim)
arguments
	numericdata
end
arguments (Repeating)
	dim			
	ind4EachDim	
end
% if dim is a scalar, ind4eachdim is a vector
% if dim is a vector, ind4eachdim is a cell of vectors
IND=repmat({':'},1,ndims(numericdata));
if ~isempty(dim)
	dim=cell2mat(dim);
	IND(dim)=ind4EachDim;
end
end