function plotHotEncoded(HotEncodedData,OrgData,edges,ntl)
arguments
	HotEncodedData
	OrgData
	edges
	ntl=[]
end
if isempty(ntl)
	ntl=nested_tiledlayout(2,1);
	ntl.nexttile();
end

h=plot(HotEncodedData,LineWidth=2.5);
change_plot_color(h,parula(size(HotEncodedData,2)));
legend("bins "+string(1:numel(edges)-1));
grid on;

ntl.nexttile();


x=[1:size(OrgData,1),NaN]';
y=OrgData;
y(end+1,:)=NaN;

discData=discretize(OrgData,edges);
cmat=parula(numel(edges)-1);
colorarray=cmat(discData,:);
colorarray(end+1,:)=[1,0,0];

patch( ...
  'XData', x, ...
  'YData', y, ...
  'ZData', zeros(size(x)), ...
  'FaceVertexCData', colorarray, ...   % n×3 truecolor per vertex
  'FaceColor','none', ...
  'EdgeColor','interp', ...
  'LineWidth',2 ...
);
hold on;
modedges=edges(~isinf(edges));
yline(modedges,LineWidth=1.5);
grid on;
hold off

ntl.change_settings(LinkAxes="x")


end