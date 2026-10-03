function [taubins,taubincenter]=makeTauBins(taubincount,taumaxdeg)
arguments
	taubincount=5
	taumaxdeg=15
end

% phimaxdeg=10;
taumax=taumaxdeg*pi/180;
% phibincount=5;

taubins=linspace(-taumax,taumax,taubincount+1);
taubincenter=diff(taubins)/2+taubins(1:end-1);

taubins(1)=-Inf;
taubins(end)=Inf;

taubins=taubins(:);
taubincenter=taubincenter(:);


end