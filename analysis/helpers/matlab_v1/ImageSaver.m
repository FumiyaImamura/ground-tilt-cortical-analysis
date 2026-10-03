classdef ImageSaver < matlab.mixin.Copyable
    %IMAGESAVER Utility for saving figures as vector SVG or FIG while preserving transparency
	%% example usage
	% 	ImageSaver(gcf(),path="path...\name of file with no extension").dsvg()
	%%

	properties (Constant)
		tempPath="/path/to/data/tempfolder"

	end
    properties
        Figure  % Handle to the MATLAB figure
        Path    % Base file path (without extension)
    end

    methods
		function obj = ImageSaver(fig, options)
            arguments
                fig    matlab.ui.Figure
                options.path   string =""
				options.tempname (1,1) string ="temp"
            end
            obj.Figure = fig;
            obj.Path   = options.path;
			if options.path==""
				obj.Path=fullfile(obj.tempPath,options.tempname);
			end
        end

		function dsvg(obj)
			arguments
				obj
			end
            %DSVG Export the figure as an SVG with true vector transparency
            %   This uses exportgraphics (R2020a+) to preserve alpha channels.

			% Configure figure properties for transparent export
			fig = obj.Figure;
			% set(fig, ...
			% 	'Renderer',        'opengl', ...   % Needed for alpha support
			% 	'Color',           'none', ...     % Transparent background
			% 	'InvertHardcopy',  options.InvertHardcopy); 
			
			print(fig,obj.Path+".svg","-dsvg");

			if false
				%% unit test
				figure;tiledlayout('flow');nexttile;
				plotMeanSE(rand(10),1);nexttile;
				plotMeanSE(rand(10),1);
				ImageSaver(gcf(),tempname="aa").dsvg()
			end
		end

		function svg_modified(obj)
			fig = obj.Figure;
			plot2svg(obj.Path+".svg",fig);
		end

		function emf(obj)
			fig = obj.Figure;
			set(fig, 'Renderer', 'painters');
			print(fig, obj.Path+".emf", '-dmeta');
		end

        function fig(obj)
            %SAVEFIG Save the figure as a MATLAB FIG file
            saveas(obj.Figure, obj.Path + ".fig");
		end

		function pdf(obj)
			fig=obj.Figure;
			% 1) Make figure background transparent & preserve on‐screen look
			set(fig, ...
				'Color',          'none', ...   % transparent figure background
				'InvertHardcopy', 'off', ...    % don’t flip colors when printing
				'Renderer',       'opengl');    % ensures alpha support

			% 2) Export to PDF as pure vector with embedded transparency
			exportgraphics(fig, obj.Path+".pdf", ...
				'ContentType',    'vector', ...
				'BackgroundColor','none');
		end
		function png(obj,options)
			arguments
				obj
				options.Resolution=300
			end
			fig=obj.Figure;
			exportgraphics(fig,obj.Path+".png",Resolution=options.Resolution);
		end
    end
end
