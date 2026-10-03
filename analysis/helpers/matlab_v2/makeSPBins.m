function [edges_cellVector,bincenter_cellVector]=makeSPBins(options)
arguments
	options.phibincount=5
	options.taubincount=6
	options.Xbincount=1
	options.Ybincount=1
	options.Zbincount=5
	options.Xrange_in_mm=[0,0]
	options.Yrange_in_mm=[0,0]
	options.Zrange_in_mm=[-5,5]
end

edges_cellVector=[];
bincenter_cellVector=[];
edges=[];
bincenter=[];

i=1;
[edges{i},bincenter{i}]=makeThetaBins();
i=2;
[edges{i},bincenter{i}]=makePhiBins(options.phibincount);
i=3;
[edges{i},bincenter{i}]=makeTauBins(options.taubincount);
i=4;
[edges{i},bincenter{i}]=makeBinsFromData(options.Xbincount,options.Xrange_in_mm);
i=5;
[edges{i},bincenter{i}]=makeBinsFromData(options.Ybincount,options.Yrange_in_mm);
i=6;
[edges{i},bincenter{i}]=makeBinsFromData(options.Zbincount,options.Zrange_in_mm);

edges_cellVector=edges;
bincenter_cellVector=bincenter;

end

function [edges,bincenter]=makeBinsFromData(bincount,data_range)
if data_range(1)>data_range(2)
	error('range_in_mm needs to be [min,max]');
end
edges=linspace(data_range(1),data_range(2),bincount+1);
bincenter=edges(1:end-1)+diff(edges,1)/2;
edges(1)=-Inf;
edges(end)=Inf;

edges=edges(:);
bincenter=bincenter(:);
end