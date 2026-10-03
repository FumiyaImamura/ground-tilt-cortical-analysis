function [phiedges,phibincenter]=makePhiBins(phibincount,phimaxdeg)
arguments
	phibincount=5
	phimaxdeg=10
end

% phimaxdeg=10;
phimax=phimaxdeg*pi/180;
% phibincount=5;

% create phibins
phiedges=[0:phimax/phibincount:phimax];
phibincenter=phiedges+phimax/phibincount/2;
phibincenter=phibincenter(1:end-1);
phiedges(end)=Inf;
phiedges=phiedges(:);
phibincenter=phibincenter(:);

end