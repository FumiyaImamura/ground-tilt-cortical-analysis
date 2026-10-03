function [mean_handle,se_handle]=plotMeanSE(data, dim, options)
% plotMeanSE plots the mean of data along the time dimension with
% a shaded region representing the standard error.
%
%   plotMeanSE(data, dim)
%   plotMeanSE(data, dim, options)
%
%   data    : A 2D numeric matrix.
%   dim     : The dimension corresponding to time.
%             (1: time along rows, 2: time along columns)
%   options : (Optional) A structure with fields:
%             - t : A vector of time values for the x-axis.
%
% Example:
%   data = randn(100, 20);
%   options.x = linspace(0, 10, 100);
%   plotMeanSE(data, 1, options);

    arguments
        data (:,:) {mustBeNumeric}
        dim (1,1) {mustBeNumeric, mustBeMember(dim, [1, 2])}
        options.x = []
		options.Color=[1,1,1]*0.5
		options.ShadeAlpha=0.5
		options.LineWidth=2
		options.DispName_Mean=[]
		options.DispName_SE=[]
    end

    % Determine the replicate dimension (the one that is not the time dimension)
    repDim = 3 - dim;  % works for a 2D matrix
    
    % Compute the mean and standard error along the replicate dimension
    m  = mean(data, repDim);
    se = std(data, 0, repDim) / sqrt(size(data, repDim));
    
    % Ensure m and se are row vectors for consistent plotting
    m  = m(:)';
    se = se(:)';
    
    % Create or validate the time vector using options.x
    if isempty(options.x)
        t = 1:length(m);
    else
        t = options.x;
        t = t(:)'; % ensure t is a row vector
        if length(t) ~= length(m)
            error('Time vector length must match the number of time points in the data.');
        end
    end
    
    % Compute the upper and lower bounds for the shaded area
    upper = m + se;
    lower = m - se;
    
    % Plot the shaded error region
	% if isempty(options.Color)
	% 	options.Color=[1,1,1]*0.5;
	% end
	% if isempty(options.ShadeAlpha)
	% 	options.ShadeAlpha=0.5;
	% end
	Log=(~isnan(t)).*(~isnan(upper)).*(~isnan(lower));
	Log=logical(Log);
	T=t(Log);
	Upper=upper(Log);
	Lower=lower(Log);

    holdlog=ishold;
    % Plot the mean as a line
	if isempty(options.DispName_Mean)
		mean_handle=plot(t, m, Color=options.Color, LineWidth=options.LineWidth,LineStyle="-");
		mean_handle.Annotation.LegendInformation.IconDisplayStyle="off";
	else
		mean_handle=plot(t, m, Color=options.Color, LineWidth=options.LineWidth,LineStyle="-",DisplayName=options.DispName_Mean);
	end
	hold on;
	if isempty(options.DispName_SE)
		se_handle=fill([T, fliplr(T)], [Upper, fliplr(Lower)], options.Color, ...
			'EdgeColor', 'none', 'FaceAlpha', options.ShadeAlpha);
		se_handle.Annotation.LegendInformation.IconDisplayStyle="off";
	else
		se_handle=fill([T, fliplr(T)], [Upper, fliplr(Lower)], options.Color, ...
			'EdgeColor', 'none', 'FaceAlpha', options.ShadeAlpha,DisplayName=options.DispName_SE);
	end
    hold off;

	modes = {'off','on'};     % modes{1} = 'off', modes{2} = 'on'
	hold( modes{ holdlog+1 } );  % ishold returns 0 or 1 → index 1 or 2

    
    % Label the plot
    % xlabel('Time');
    % ylabel('Mean');
    % title('Mean and Standard Error');
    grid on;
end
