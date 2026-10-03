function change_plot_color(plot_handle,color_array,alpha)
arguments
	plot_handle
	color_array
	alpha=1;
end
color_array(:,end+1)=alpha;
set(plot_handle, {'Color'}, mat2cell(color_array, ones(1, size(color_array,1)), size(color_array,2)))
end