function polarTicks8direction()
% this function assumes that the angle 0 corresponds to datapoint between
% direction 1 (45 degrees) and direction 8 (360 degrees)
deg_list=45/2+45*[0:7];
thetaticks(deg_list);
thetaticklabels(compose("%d",[1:8]*45)+"\circ")

end