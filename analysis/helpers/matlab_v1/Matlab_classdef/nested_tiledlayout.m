classdef nested_tiledlayout < handle
	properties
		handles = {}
		tile_indices = {}
		all_handles_layer = {}
		all_axes=gobjects(0)
	end
	methods
		%%
		% Constructor: Initializes the nested tiled layout object and
		% creates a figure and the first layer of tiledlayout
		% If row and col are not specified, a 'flow' layout is used.
		function obj = nested_tiledlayout(row, col, options)
			arguments
				row = []
				col = []
				options.?matlab.graphics.layout.TiledChartLayout
			end
			if ~isfield(options, "Padding")
				options.Padding = "tight";
			end
			if ~isfield(options, "TileSpacing")
				options.TileSpacing = "tight";
			end
			nvCell = namedargs2cell(options);
			Fig = figure;
			Fig.Units = "normalized";
			% Fig.Position = [0, 0, 1, 1];
			Fig.WindowState="maximized";
			try
				t = tiledlayout(row, col, nvCell{:});
			catch
				t = tiledlayout('flow', nvCell{:});
			end
			obj.handles{1} = t;
			obj.tile_indices{1} = 0;
			obj.all_handles_layer{1} = {t};
		end
		%%
		% make_tiles: Creates a new nested tiled layout within the current layout.
		% When row and col are not specified, it uses a 'flow' layout.
		function make_tiles(obj, row, col, options)
			arguments
				obj
				row = []
				col = []
				options.?matlab.graphics.layout.TiledChartLayout
			end
			if ~isfield(options, "Padding")
				options.Padding = "tight";
			end
			if ~isfield(options, "TileSpacing")
				options.TileSpacing = "tight";
			end
			nvCell = namedargs2cell(options);
			obj.tile_indices{end} = obj.tile_indices{end} + 1;
			if isempty(obj.handles)
				try
					t = tiledlayout(row, col, nvCell{:});
				catch
					t = tiledlayout('flow', nvCell{:});
				end
			else
				try
					t = tiledlayout(obj.handles{end}, row, col, nvCell{:});
				catch
					t = tiledlayout(obj.handles{end}, 'flow', nvCell{:});
				end
				t.Layout.Tile = obj.tile_indices{end};
			end
			obj.handles{end+1} = t;
			obj.tile_indices{end+1} = 0;
			if length(obj.all_handles_layer) < length(obj.handles)
				obj.all_handles_layer{length(obj.handles)} = {};
			end
			obj.all_handles_layer{length(obj.handles)}{end+1} = t;
		end
		%%
		% nexttile: Advances to the next tile in the current tiled layout.
		function ax=nexttile(obj,tilelocation,span)
			arguments
				obj
				tilelocation=[]
				span=[]
			end
			if isempty(tilelocation)
				if isempty(span)
					ax= nexttile(obj.handles{end});
				else
					ax= nexttile(obj.handles{end},span);
				end
			else
				if isempty(span)
					ax= nexttile(obj.handles{end},tilelocation);
				else
					ax= nexttile(obj.handles{end},tilelocation,span);
				end
			end
			obj.all_axes(end+1)=ax;
			obj.tile_indices{end} = obj.tile_indices{end} + 1;
		end
		%%
		% out_tiles: Exits the current nested tiled layout, returning to the previous layout.
		function out_tiles(obj)
			obj.handles = obj.handles(1:length(obj.handles)-1);
			obj.tile_indices = obj.tile_indices(1:length(obj.tile_indices)-1);
		end
	end
	properties
		linked_axes = {}
	end
	methods
		%%
		% change_settings: Updates various properties (titles, labels, legends, axis linking) of the current tiled layout.
		function change_settings(obj, options)
			arguments
				obj
				options.handle = []
				options.Title = []
				options.TitleSize = []
				options.childrenTitle = []		% For all axes and tiledlayout contained inside the current layer of tiledlayout
				options.childrenTitleSize = []	% For all axes and tiledlayout contained inside the current layer of tiledlayout
				options.Xlabel = []
				options.XlabelSize = []
				options.Ylabel = []
				options.YlabelSize = []
				options.Legend = []
				options.LegendSize = []
				options.LegendLoc {mustBeMember(options.LegendLoc,{'','north','east','west','south'})} = ''
				options.LinkAxes {mustBeMember(options.LinkAxes, {'xyz','x','y','z','xy','xz','yz','off'})} = 'off'
				options.ColorbarLoc {mustBeMember(options.ColorbarLoc,{'','north','east','west','south'})} = ''
				options.Prop2Change	={}% 1d cell array. first element is a string specifying the property name. Set this to string array if it is a nested property.
									% the second element is the value you
									% want to set.
									% keep alternating the [property names]
									% and [property value to set]
				options.childrenProp2Change={}
				options.LinkAxes_all {mustBeMember(options.LinkAxes_all, {'xyz','x','y','z','xy','xz','yz','off'})} = 'off'
				options.CLim (1,:) =[];
			end
			if ~isempty(options.handle)
				tl = options.handle;
			else
				tl = obj.handles{end};
			end
			if ~isempty(options.Title)
				tl.Title.String = options.Title;
			end
			if ~isempty(options.TitleSize)
				tl.Title.FontSize = options.TitleSize;
			end
			if ~isempty(options.childrenTitle)
				List = find_specified_children(tl, "axes");
				if ~isempty(List)
					arrayfun(@(a) title(a, options.childrenTitle), List);
				end
				List = find_specified_children(tl,"tiledlayout");
				if ~isempty(List)
					arrayfun(@(a) title(a, options.childrenTitle), List);
				end
			end
			if ~isempty(options.childrenTitleSize)
				List = find_specified_children(tl, "axes");
				if ~isempty(List)
					arrayfun(@(a) set(get(a, "Title"), "FontSize", options.childrenTitleSize), List);
				end
				List = find_specified_children(tl,"tiledlayout");
				if ~isempty(List)
					arrayfun(@(a) set(get(a, "Title"), "FontSize", options.childrenTitleSize), List);
				end
			end
			if ~isempty(options.Xlabel)
				tl.XLabel.String = options.Xlabel;
			end
			if ~isempty(options.XlabelSize)
				tl.XLabel.FontSize = options.XlabelSize;
			end
			if ~isempty(options.Ylabel)
				tl.YLabel.String = options.Ylabel;
			end
			if ~isempty(options.YlabelSize)
				tl.YLabel.FontSize = options.YlabelSize;
			end
			if ~isempty(options.Legend)
				List = find_specified_children(tl, "axes");
				arrayfun(@(a) legend(a, options.Legend), List);
			end
			if ~isempty(options.LegendSize)
				List = find_specified_children(tl, "legend");
				arrayfun(@(a) set_value(a, "FontSize", options.LegendSize), List);
			end
			if ~isempty(options.LegendLoc)
				List = find_specified_children(tl, "legend");
				strlist = arrayfun(@(a) a.String, List, 'UniformOutput', false);
				% if ~isequal(strlist{:}); error("Legends must be all the same"); end
				List(1).Layout.Tile = options.LegendLoc;
				arrayfun(@(a) set_value(a, "Visible", "off"), List(2:end));
			end
			if ~isempty(options.ColorbarLoc)
				List = find_specified_children(tl, "colorbar");
				List(1).Layout.Tile = options.ColorbarLoc;
				arrayfun(@(a) set_value(a, "Visible", "off"), List(2:end));
			end
			if string(options.LinkAxes) ~= "off"
				List = find_specified_children(tl, "axes");
				if string(options.LinkAxes) == "xyz"
					obj.linked_axes{end+1} = linkprop(List, {'CameraPosition', 'CameraTarget', 'CameraUpVector'});
				end
				linkaxes(List, options.LinkAxes);
			end
			if ~isempty(string(options.CLim))
				List = find_specified_children(tl, "axes");
				for i=1:length(List)
					set(List(i),"CLim",options.CLim);
				end
			end
			if ~isempty(options.Prop2Change)
				obj.setGraphicProperties(tl,options.Prop2Change);
			end
			if ~isempty(options.childrenProp2Change)
				currhandle=tl.Children;
				obj.setGraphicProperties(currhandle,options.childrenProp2Change);
			end
			if string(options.LinkAxes_all)~="off"
				List=findobj(tl,"Type","axes");
				List = List(~ismember(get(List,'Tag'), {'legend','Colorbar'}));
				if string(options.LinkAxes_all) == "xyz"
					obj.linked_axes{end+1} = linkprop(List, {'CameraPosition', 'CameraTarget', 'CameraUpVector'});
				end
				linkaxes(List, options.LinkAxes_all);
			end
		end
		%%
		% change_settings_layer: Applies changes to settings for all tiled layouts in a specified layer.
		% This can be executed retrospectively
		function change_settings_layer(obj, layer2change, options)
			arguments
				obj
				layer2change
				options.Title = []
				options.TitleSize = []
				options.childrenTitle = []
				options.childrenTitleSize = []
				options.Xlabel = []
				options.XlabelSize = []
				options.Ylabel = []
				options.YlabelSize = []
				options.Legend = []
				options.LegendSize = []
				options.LegendLoc {mustBeMember(options.LegendLoc,{'','north','east','west','south'})} = ''
				options.LinkAxes {mustBeMember(options.LinkAxes, {'xyz','x','y','z','xy','xz','yz','off'})} = 'off'
				options.ColorbarLoc {mustBeMember(options.ColorbarLoc,{'','north','east','west','south'})} = ''
				options.Prop2Change	={}% 1d cell array. first element is a string specifying the property name. Set this to string array if it is a nested property.
									% the second element is the value you
									% want to set.
									% keep alternating the [property names]
									% and [property value to set]
				options.childrenProp2Change={}
				options.LinkAxes_all {mustBeMember(options.LinkAxes_all, {'xyz','x','y','z','xy','xz','yz','off'})} = 'off'
			end
			handlelist = obj.all_handles_layer{layer2change};
			for i = 1:length(handlelist)
				opts = options;
				opts.handle = handlelist{i};
				opts = namedargs2cell(opts);
				obj.change_settings(opts{:});
			end
		end
		%%
	end
	methods (Access=private)
		function setGraphicProperties(obj, h, Prop2Change)
			% SETGRAPHICPROPERTIES Set properties on a graphics object (or array of objects),
			% including nested properties.
			%
			%   setGraphicProperties(obj, h, Prop2Change) sets the properties of the graphics
			%   object(s) h as specified in Prop2Change. This cell array should alternate between
			%   property name specifiers and property values. A property name specifier is:
			%       - a char or scalar string for a simple property, or
			%       - a string array for a nested property (e.g. ["Title","String"]).
			%
			%   Example:
			%       % Set the 'CLim' property and a nested 'Title.String' property:
			%       Prop2Change = { 'CLim', [0 1], ["Title","String"], 'New Title' };
			%       setGraphicProperties(obj, h, Prop2Change);
			%
			%   Note: When a property (or nested property) refers to an array of graphics objects,
			%   the function automatically iterates over each element.

			props = Prop2Change;

			% Ensure there is an even number of elements.
			if mod(length(props), 2) ~= 0
				error('Prop2Change must have an even number of elements: alternating property names and values.');
			end

			% Loop through each object in the handle array.
			for k = 1:numel(h)
				curObj = h(k);

				% Loop over each property-value pair.
				for i = 1:2:length(props)
					propNameSpec = props{i};
					newValue = props{i+1};

					if ischar(propNameSpec) || (isstring(propNameSpec) && numel(propNameSpec)==1)
						% Simple (non-nested) property:
						try
							% Use the built-in set function
							set(curObj, propNameSpec, newValue);
						catch ME
							warning('Failed to set property "%s" on object %d: %s', char(propNameSpec), k, ME.message);
						end

					elseif isstring(propNameSpec) && numel(propNameSpec)>1
						% Nested property chain. Use recursion so that if any step of the chain
						% returns an array of objects, we iterate through each element.
						try
							setNestedProperty(curObj, propNameSpec, newValue);
						catch ME
							warning('Failed to set nested property chain "%s" on object %d: %s', ...
								strjoin(propNameSpec, '.'), k, ME.message);
						end

					else
						error('Invalid property name specifier. Use a char for a simple property or a string array for nested properties.');
					end
				end
			end
			function setNestedProperty(objOrArray, propChain, newValue)
				% setNestedProperty Recursively traverse and set a nested property on an object or array.
				%
				%   setNestedProperty(objOrArray, propChain, newValue) takes an object (or an array of objects)
				%   'objOrArray' and a string array 'propChain' that specifies the nested property path.
				%   It iterates through every element in any array of objects it encounters.
				%
				%   The recursion stops when propChain has only one element; in that case the property is
				%   set directly on all objects in objOrArray.

				% Make sure we are working with an array (even if a scalar, numel returns 1)
				for idx = 1:numel(objOrArray)
					currentObj = objOrArray(idx);
					if numel(propChain)==1
						% Base case: set the property on currentObj.
						try
							set(currentObj, propChain(1), newValue);
						catch ME
							warning('Failed to set nested property "%s": %s', char(propChain(1)), ME.message);
						end
					else
						% Recursive case: Get the next part of the chain.
						try
							nextObj = currentObj.(propChain(1));
						catch ME
							warning('Failed to access nested property "%s": %s', char(propChain(1)), ME.message);
							continue;
						end
						% Recursively call on the remainder of the chain.
						setNestedProperty(nextObj, propChain(2:end), newValue);
					end
				end
			end
		end
	end
end

%%
% find_specified_children: Returns a list of children of a given type from a tiled layout.
function newlist = find_specified_children(tl_handle, graphic_type)
List = tl_handle.Children;
newlist = gobjects();
counter = 1;
for i = 1:length(List)
	if string(List(i).Type) == string(graphic_type)
		newlist(counter) = List(i);
		counter = counter + 1;
	end
end
if counter==1
	newlist=[];
end
end

% set_value: Sets a specified property value for a given graphics object.
function variable = set_value(variable, prop, value)
variable.(string(prop)) = value;
end
