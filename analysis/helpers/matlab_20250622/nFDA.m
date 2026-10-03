function [coeff,Score,A,B,ScoreA,ScoreB,val]=nFDA(neuw,arg,w)
% neuw * coeff = Score

% nFDA - Perform dimensionality reduction on neural data using various methods.
%
% Syntax: [coeff, Score, A, B, ScoreA, ScoreB, val] = nFDA(neuw, arg, w)
%
% Inputs:
%   neuw - A 3D matrix representing neural data with dimensions(time within trial * cells * trials)
%         =(t, nCells, nTrials).
%   arg  - A string specifying the method for dimensionality reduction. 
%          Options are "nFDA", "Rayleigh", or "TRCA".
%   w    - A scalar weight parameter (0 to 1) used only when arg is "nFDA". 
%          It balances the trade-off between maximizing within-trial variance 
%          and minimizing across-trial variance. Smaller values of w prioritize 
%          within-trial variance, while larger values prioritize across-trial variance.
%          Around 0.5~0.85 is a good range from experience
%
% Outputs:
%   coeff  - A 2D matrix of coefficients for dimensional reduction. The first 
%            column contains coefficients for the first axis.
%   Score  - A 3D matrix with the same dimensions as neuw, containing the scores 
%            for each axis, similar to the output of the PCA function in MATLAB.
%   A      - Matrix storing the results of within-trial variance for each axis.
%   B      - Matrix storing the results of across-trial variance for each axis.
%   ScoreA - A vector containing within-trial variability scores for each axis.
%   ScoreB - A vector containing across-trial variability scores for each axis.
%   val    - Eigenvalues for each methods.
%

%
% Methods:
%   - "nFDA": Performs dimension reduction with the specified weight w.
%   - "Rayleigh": Maximizes the ratio of within-trial variance to across-trial variance.
%   - "TRCA": Another dimension reduction method for comparison. Confer:https://pubmed.ncbi.nlm.nih.gov/22922468/

nTrial=size(neuw,3);
nTime=size(neuw,1);
% make covariance matrix for within trial variance
IavgA=eye(nTime,nTime)-(1/nTime)*ones(nTime,nTime);
IxA=pagemtimes(IavgA,neuw);
A=pagemtimes(permute(neuw,[2 1 3]),IxA);
Araw=sum(A,3);
Araw=(1/2)*(Araw+Araw'); % make it hermitian (theoretically it is hermitian, but it is not because of the floating point error)
% plus hermitian has much faster algorithm for calculating eigenvalues
A=(1/(nTrial*nTime))*Araw;
A=A./trace(A);


% make covariance matrix for across trial variance
IavgB=eye(nTrial,nTrial)-(1/nTrial)*ones(nTrial,nTrial);
IxB=pagemtimes(IavgB,permute(neuw,[3 2 1]));

%%%% 2024/09/24 %%%%
% the output array was too big

%%%%%%%%% previous code
% B=pagemtimes(permute(neuw,[2 3 1]),IxB);
% Braw=sum(B,3);
%%%%%%%%%

%%%%%%%%% new code
pneuw=permute(neuw,[2,3,1]);
temp=pneuw(:,:,1)*IxB(:,:,1);
Braw=zeros(size(temp));

for i=1:size(IxB,3)
	temp=pneuw(:,:,i)*IxB(:,:,i);
	Braw=Braw+temp;
	% disp(append(num2str(i),' / ',num2str(size(IxB,3))));
end

%%%%%%%%%

Braw=(1/2)*(Braw+Braw');
B=(1/(nTrial*nTime))*Braw;
B=B./trace(B);

% solve by difference
if strcmpi(arg,'nFDA')

	switch nargin
		case 3
			w1=w; % parameter to change the weight of A or B
			[eigv,val]=eig((1-w1)*A-w1*B);
			[val,ind]=sort(diag(real(val)),1,'descend');
			eigv=eigv(:,ind);
			eigv=eigv./(sqrt(dot(eigv,eigv,1)));
			%get scores of each axis
			Score=pagemtimes(neuw,eigv);

			ScoreA=dot(eigv,A*eigv,1);
			ScoreB=dot(eigv,B*eigv,1);
			coeff=eigv;
		otherwise
			error('Not enough inputs. Specify parameter ''w''.')
	end


elseif strcmpi(arg,'Rayleigh')
	switch nargin
		case 2
			% rayleigh quotient
			[eigv2,val2]=eig(A,B);
			[val2,ind2]=sort(diag(val2),'descend');
			eigv2=eigv2(:,ind2);
			eigv2=eigv2./(sqrt(dot(eigv2,eigv2,1)));

			%get scores of each axis
			Score=pagemtimes(neuw,eigv2);

			ScoreA=real(dot(eigv2,A*eigv2,1))';
			ScoreB=real(dot(eigv2,B*eigv2,1))';

			coeff=eigv2;
			val=val2;
		otherwise
			error('Too many inputs')
	end


elseif strcmpi(arg,'TRCA')
	switch nargin
		case 2
			%TRCA
			neuwsum=sum(neuw,3);
			Craw=neuwsum'*IavgA*neuwsum;
			Craw=(1/2)*(Craw+Craw');
			C=Craw-A;
			Cnorm=Craw/((nTrial^2)*nTime);

			[eigv3,val3]=eig(C,A);
			[val3,ind3]=sort(diag(val3),'descend');
			eigv3=eigv3(:,ind3);
			%normalize coefficients
			eigv3=eigv3./(sqrt(dot(eigv3,eigv3,1)));
			%get scores of each axis
			Score=pagemtimes(neuw,eigv3);

			ScoreA=real(dot(eigv3,A*eigv3,1));
			ScoreB=real(dot(eigv3,B*eigv3,1));
			coeff=eigv3;
			val=val3;
		otherwise
			error('Too many inputs')
	end
else
	error('Not enough inputs or incorrect input. 2nd input should be ''nFDA'', ''rayleigh'', or ''TRCA''.')
end
coeff=real(coeff);
Score=real(Score);
end