function c=FastConv(a,b,dim,options)
arguments
	a
	b
	dim
	options.shape {mustBeMember(options.shape,["full","same_as_a"])}="same_as_a"
	options.normalize (1,1) logical =false
	options.mirror_endpoints (1,1) logical=false
end

if options.mirror_endpoints
	length2flip=min(size(a,dim),floor(size(b,dim)/2));
	IND=repmat({':'},1,max(ndims(a),dim));
	IND{dim}=1:length2flip;
	pre_a=a(IND{:});
	fpre_a=flip(pre_a,dim);
	IND{dim}=size(a,dim)-length2flip+1:size(a,dim);
	post_a=a(IND{:});
	fpost_a=flip(post_a,dim);
	a=cat(dim,fpre_a,a,fpost_a);
end


% Calculates normalized crosscorrelation of a and b along dim
% "a" and "b" can be any vector or high dimensional matrix but there sizes must
% be compatible except for the "dim" dimension. Compatible meaning the
% sizes of "a" and "b" along all dimensions are either the same or one of them is 1
% The output "c" is the same size as a.*b except for dim dimension which has the size of size(a,dim)+size(b,dim)-1.
% The first value along the dim dimension represents convolution between
% first elements of "a" along dim dimension and last elements of "b" along dim
% dimension.
sa=size(a);maxdima=length(sa);dimveca=(1:maxdima);dimveca(dim)=[];
sb=size(b);maxdimb=length(sb);dimvecb=(1:maxdimb);dimvecb(dim)=[];
ap=permute(a,[dim,dimveca]);sap=size(ap);
bp=permute(b,[dim,dimvecb]);sbp=size(bp);
anew=cat(1,zeros([size(bp,1)-1,sap(2:end)]),ap);
bnew=cat(1,bp,zeros([size(ap,1)-1,sbp(2:end)]));
c1=ifft(fft(anew,[],1).*conj(fft(bnew,[],1)),[],1);
if options.normalize
	d1=sqrt(movsum(anew.^2,[0,size(bp,1)-1],1));
	d2=sqrt(movsum(flip(bnew,1).^2,[0,size(ap,1)-1],1));
	c1=c1./(d1.*d2);
end
M=max(maxdima,maxdimb);
ind=1:M;
ind1=circshift(1:dim,-1);
ind(1:dim)=ind1;
c=permute(c1,ind);
if options.shape=="same_as_a"
	% Extract the central part along dimension "dim" so that c has the same size as a:
	startIdx = floor((size(b,dim)-1)/2) + 1;
	idx = repmat({':'}, 1, ndims(c));
	idx{dim} = startIdx:(startIdx + size(a,dim) - 1);
	c = c(idx{:});
end

if options.mirror_endpoints
	IND=repmat({':'},1,max(ndims(a),dim));
	IND{dim}=length2flip+1:size(a,dim)-length2flip;
	c=c(IND{:});
end


end
