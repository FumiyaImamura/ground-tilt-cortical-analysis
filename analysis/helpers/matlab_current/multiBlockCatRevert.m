function x=multiBlockCatRevert(newx,info)

dim2concat=info.dim2concat;
dim2squish_vec=info.dim2squish_vec;
dim2squish_vec=dim2squish_vec(end:-1:1);	% reverse the order
sizeEachStep=info.sizeEachStep;
sizeEachStep=sizeEachStep(end:-1:1,:);		% reverse the order

x=newx;
for i=1:length(dim2squish_vec)
	x=BlockCatRevert(x,dim2squish_vec(i),dim2concat,sizeEachStep(i,:));
end

if false
	%% test code
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