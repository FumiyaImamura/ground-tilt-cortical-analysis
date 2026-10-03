function xmat = plotSideBySide(data, options)
% plotSideBySide Plots multiple slices of data side-by-side with custom options.
%
%   xmat = plotSideBySide(data, namestr, options) plots each 2D slice of the
%   3D numeric array 'data' side-by-side with x-axis offsets. It uses the
%   provided 'namestr' for labeling each slice and 'options' for custom plot
%   properties. If options.xmat is provided, it is used as the x-axis values.
%
%   Parameters:
%     data    - A 3D numeric array of size (m x n x p) where:
%               m: number of data points per curve,
%               n: number of curves per slice,
%               p: number of slices.
%
%     namestr - A string array of size (p x 1) containing names for each slice.
%
%     options - A structure with the following fields:
%               Xlabel    : Text label for the x-axis.
%               Ylabel    : Text label for the y-axis.
%               Title     : Plot title.
%               Color     : Line color (default: 'red').
%               Linewidth : Width of the plot lines (default: 1).
%               xmat      : Custom x-axis values for each slice (default: []).
%
%   Returns:
%     xmat    - A matrix containing the x-axis values used for each slice.
arguments
	data (:,:,:) {mustBeNumeric}
	options.xticklabels (:,1) string=[]
	options.Xlabel string=[]
	options.Ylabel string=[]
	options.Title string=[]
	options.Color =[1,1,1]*0.6
	options.Linewidth =1
	options.xmat (:,:) {mustBeNumeric} =[];
	options.PlotMeanSE (1,1) logical =true
	options.ShadeAlpha=0.5
	options.xtickFontSize=[]
	options.plotIndividualData=true;
	options.MeanColor="blue"
end
% Parse the options structure and assign variables.
Xlabel    = options.Xlabel;
Ylabel    = options.Ylabel;
Title     = options.Title;
Color     = options.Color;
Linewidth = options.Linewidth;
xmat      = options.xmat;

% Example input usage:
% data: A 3D array of size (m x n x p)
% xticklabels: A string array of size (1 x p)

% Get the size of the input data.
[m, n, p] = size(data);

% Initialize the x-axis offset.
x_offset = 0;

washold=ishold(gca());

% Open a new figure and prepare to hold multiple plots.

% Check if custom x-axis matrix is provided.
if isempty(xmat)
    % If no custom x-axis matrix, create one for each slice.
    xmat = [];
    for k = 1:p
        % Calculate x-axis values for the current slice with an offset.
        x = (1:m) + x_offset;  % Create a shifted x-axis for this slice.
        xmat(:, k) = x;        % Store the x-axis values in the xmat matrix.

        % Plot every curve (each column) in the current slice.
        % for j = 1:n
        %     plot(x, data(:, j, k), 'Color', Color, 'LineWidth', Linewidth);
		% end
		if options.plotIndividualData
			plot(x, data(:, :, k), 'Color', Color, 'LineWidth', Linewidth);
			hold on;
		end
		if options.PlotMeanSE
			plotMeanSE(data(:, :, k),1,Color=options.MeanColor,ShadeAlpha=options.ShadeAlpha,x=x);
			hold on;
		end

        % Update the x-axis offset to separate this slice from the next.
        x_offset = x_offset + m + 1;  % Increment offset by the number of data points plus a gap.
    end
else
    % If a custom x-axis matrix is provided, use it.
    for k = 1:p
        % Use the provided x-axis values for this slice.
        x = xmat(:, k);

        % Plot every curve (each column) in the current slice.
        for j = 1:n
            plot(x, data(:, j, k), 'Color', Color, 'LineWidth', Linewidth);
			hold on;
        end

        % Update the x-axis offset (still calculated for overall plot limits).
        x_offset = x_offset + m + 1;
    end
end

% Set the x-axis ticks at the midpoint of each slice.
x_ticks = m / 2 + (0:p-1) * (m + 1);
xticks(x_ticks);


% Label each tick with the corresponding name from xticklabels.
if ~isempty(options.xticklabels)
xticklabels(options.xticklabels);

if ~isempty(options.xtickFontSize)
	ax=gca;
	ax.XAxis.FontSize=options.xtickFontSize;

% Set the x-axis label if provided.
if ~isempty(Xlabel)
    xlabel(Xlabel);
end

% Set the y-axis label if provided.
if ~isempty(Ylabel)
    ylabel(Ylabel);
end

% Set the plot title if provided.
if ~isempty(Title)
    title(Title);
end

% Optional: Uncomment the next line to add grid lines to the plot.
% grid on;

% Set the x-axis limits to span from 0 to the final offset value.
xlim([0, x_offset]);

% Release the hold on the current figure.
if ~washold
	hold off;
end

end
end
