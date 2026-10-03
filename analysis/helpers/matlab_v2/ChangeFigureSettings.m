function ChangeFigureSettings(options)
arguments
	options.MasterTitle = []      % Overall (master) title string
	options.Title = []            % Axes title string (applied per axes)
	options.Legend = []           % Legend labels (cell or string array)
	options.LegendLoc = []        % Legend location (for axes or tiled layout)
	options.MasterTitleSize = []  % Font size for master title
	options.TitleSize = []        % Font size for axes titles
	options.XLabelSize = []       % Font size for x-axis labels
	options.YlabelSize = []       % Font size for y-axis labels
	options.LegendSize = []       % Font size for legend text
	options.Linkaxes {mustBeMember(options.Linkaxes, {'xyz','x','y','z','xy','xz','yz','off'})} = 'off'
	options.Xlabel = []           % NEW: x-axis label text
	options.Ylabel = []           % NEW: y-axis label text
end

% Get the current figure
fig = gcf;

%% Locate all axes (excluding legends and colorbars)
allAxes = findobj(fig, 'Type', 'axes');
ax = allAxes(~arrayfun(@(h) any(strcmp(h.Tag,{'legend','Colorbar'})), allAxes));
numAxes = numel(ax);

%% MASTER TITLE HANDLING
% Check for tiled layout objects first.
tiledObjs = findall(fig, 'Type', 'tiledlayout');
if ~isempty(tiledObjs)
	% Use the first tiled layout as the master container.
	t = tiledObjs(1);
	if ~isempty(options.MasterTitle)
		if isprop(t, 'Title') && ~isempty(t.Title) && isgraphics(t.Title)
			t.Title.String = options.MasterTitle;
		else
			t.Title = text(t, 0.5, 1.05, options.MasterTitle, ...
				'Units', 'normalized', 'HorizontalAlignment', 'center');
		end
	end
	if ~isempty(options.MasterTitleSize)
		if isprop(t, 'Title') && ~isempty(t.Title) && isgraphics(t.Title)
			t.Title.FontSize = options.MasterTitleSize;
		end
	end
else
	% If no tiled layout exists, check for an sgtitle object.
	masterTitleObj = findall(fig, 'Type', 'text', 'Tag', 'sgtitle');
	if ~isempty(masterTitleObj)
		if ~isempty(options.MasterTitle)
			masterTitleObj.String = options.MasterTitle;
		end
		if ~isempty(options.MasterTitleSize)
			masterTitleObj.FontSize = options.MasterTitleSize;
		end
	else
		if ~isempty(options.MasterTitle)
			sg = sgtitle(options.MasterTitle);
			if ~isempty(options.MasterTitleSize)
				sg.FontSize = options.MasterTitleSize;
			end
		end
	end
end

%% INDIVIDUAL AXES TITLES AND LABELS
if numAxes == 1
	% SINGLE AXES: apply changes to that one axes.
	currAx = ax;
	axTitle = get(currAx, 'Title');
	if ~isempty(options.Title)
		axTitle.String = options.Title;
	end
	if ~isempty(options.TitleSize) && isgraphics(axTitle)
		axTitle.FontSize = options.TitleSize;
	end

	% NEW: Update XLabel and YLabel text if provided.
	if ~isempty(options.Xlabel)
		currAx.XLabel.String = options.Xlabel;
	end
	if ~isempty(options.Ylabel)
		currAx.YLabel.String = options.Ylabel;
	end

	if ~isempty(options.XLabelSize) && isgraphics(currAx.XLabel)
		currAx.XLabel.FontSize = options.XLabelSize;
	end
	if ~isempty(options.YlabelSize) && isgraphics(currAx.YLabel)
		currAx.YLabel.FontSize = options.YlabelSize;
	end
else
	% MULTIPLE AXES: update each axes individually.
	for i = 1:numAxes
		axTitle = get(ax(i), 'Title');
		if ~isempty(options.Title)
			axTitle.String = options.Title;
		end
		if ~isempty(options.TitleSize) && isgraphics(axTitle)
			axTitle.FontSize = options.TitleSize;
		end

		% NEW: Update axis label text if provided.
		if ~isempty(options.Xlabel)
			ax(i).XLabel.String = options.Xlabel;
		end
		if ~isempty(options.Ylabel)
			ax(i).YLabel.String = options.Ylabel;
		end

		if ~isempty(options.XLabelSize) && isgraphics(ax(i).XLabel)
			ax(i).XLabel.FontSize = options.XLabelSize;
		end
		if ~isempty(options.YlabelSize) && isgraphics(ax(i).YLabel)
			ax(i).YLabel.FontSize = options.YlabelSize;
		end
	end
end

%% LEGEND HANDLING
legendObjs = findobj(fig, 'Type', 'legend');
if isempty(legendObjs)
	if ~isempty(options.Legend)
		% Check if the first axes belongs to a tiled layout.
		if ~isempty(ax) && isgraphics(ax(1)) && strcmp(ax(1).Parent.Type, 'tiledlayout')
			if ~isempty(options.LegendLoc)
				% Create a legend with 'Location' set to 'none', then set its tile.
				lgd = legend(ax(1), options.Legend, 'Location', 'none');
				lgd.Layout.Tile = options.LegendLoc;
			else
				lgd = legend(ax(1), options.Legend);
			end
		else
			if ~isempty(options.LegendLoc)
				lgd = legend(ax(1), options.Legend, 'Location', options.LegendLoc);
			else
				lgd = legend(ax(1), options.Legend);
			end
		end
		if ~isempty(options.LegendSize)
			lgd.FontSize = options.LegendSize;
		end
	end
else
	for i = 1:length(legendObjs)
		if ~isempty(options.Legend)
			legendObjs(i).String = options.Legend;
		end
		if ~isempty(options.LegendSize)
			legendObjs(i).FontSize = options.LegendSize;
		end
		% For legends attached to a tiled layout, update Layout.Tile only if LegendLoc is provided.
		if ~isempty(ax) && isgraphics(ax(1)) && strcmp(ax(1).Parent.Type, 'tiledlayout')
			if ~isempty(options.LegendLoc)
				legendObjs(i).Location = 'none';
				legendObjs(i).Layout.Tile = options.LegendLoc;
			end
		else
			if ~isempty(options.LegendLoc)
				legendObjs(i).Location = options.LegendLoc;
			end
		end
	end
end
%% LINK AXES HANDLING (ADDED)
% Link axes if there are multiple axes and Linkaxes option is not 'off'
if numAxes > 1 && ~strcmp(options.Linkaxes, 'off')
	linkaxes(ax, options.Linkaxes);
end

end
