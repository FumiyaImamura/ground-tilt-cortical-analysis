function T=SPparams2TransMat(sp_array)
arguments
	sp_array
end
if size(sp_array,2)~=6
	error('sp_array needs to have 6 parameters')
end


nDims=ndims(sp_array);
IND=repmat({':'},1,nDims);

i=1;
slicer=IND;slicer{2}=i;
a=sp_array(slicer{:});
theta=a;

i=2;
slicer=IND;slicer{2}=i;
a=sp_array(slicer{:});
phi=a;

i=3;
slicer=IND;slicer{2}=i;
a=sp_array(slicer{:});
tau=a;

i=4:6;
slicer=IND;slicer{2}=i;
a=sp_array(slicer{:});
p=a;


zaxis=[sin(phi).*cos(theta),sin(phi).*sin(theta),cos(phi)];
prex=zeros(size(zaxis));
slicer=IND;slicer{2}=1;	prex(slicer{:})=cos(phi);
slicer=IND;slicer{2}=3;	prex(slicer{:})=-cos(theta).*sin(phi);
xaxis=prex./sqrt(dot(prex,prex,2));

yaxis=cross(zaxis,xaxis,2);

xaxis=cos(tau).*xaxis+sin(tau).*yaxis;
yaxis=-sin(tau).*xaxis+cos(tau).*yaxis;

R=[specialPermute(xaxis),specialPermute(yaxis),specialPermute(zaxis)];
p=specialPermute(p);

T=[R,p];
nDims=ndims(T);
IND=repmat({':'},1,nDims);
IND{1}=size(T,1)+1;
T(IND{:})=0;
IND{2}=size(T,2);
T(IND{:})=1;


end

%%
function B = specialPermute(A)
    % specialPermute rearranges the dimensions of A so that:
    %   B's 1st dim = A's 2nd dim
    %   B's 2nd dim = 1 (singleton)
    %   B's 3rd dim = A's 1st dim
    %   B's 4th dim = A's 3rd dim
    %   B's 5th dim = A's 4th dim
    %   and so on...
    
    % Number of dimensions in A
    nA = ndims(A);
    
    % Step 1: Permute dimensions so that dim 2 is first, dim 1 is second,
    % and the rest follow in ascending order.
    permOrder = [2, 1, 3:nA];
    A_temp = permute(A, permOrder);
    
    % Step 2: Insert a singleton dimension as the second dimension.
    % After permute: size(A_temp) = [size(A,2), size(A,1), size(A,3), size(A,4), ...]
    szTemp = size(A_temp);
    newSize = [szTemp(1), 1, szTemp(2:end)];
    
    % Reshape A_temp to introduce the extra singleton dimension in position 2.
    B = reshape(A_temp, newSize);
end
