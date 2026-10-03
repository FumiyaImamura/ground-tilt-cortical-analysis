function [tehtaedges,thetabincenter]=makeThetaBins
% create thetabins
tehtaedges=pi/8+2*pi*[0:8]./8;
thetabincenter=tehtaedges+pi/8;
thetabincenter=thetabincenter(1:end-1);
tehtaedges=tehtaedges(:);
thetabincenter=thetabincenter(:);
end