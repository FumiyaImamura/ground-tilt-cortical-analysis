function [mean_handle,se_handle]=plotMeanSE(data, dim, options)
% plotMeanSE plots the mean of data along the time dimension with
% either a shaded region or error bars representing the standard error.
%
%   plotMeanSE(data, dim)
%   plotMeanSE(data, dim, options)
%
%   data    : A 2D numeric matrix.
%   dim     : The dimension corresponding to time.
%             (1: time along rows, 2: time along columns)
%   options : (Optional) A structure with fields:
%             - x : A vector of time values for the x-axis.
%             - ErrorBar : true → plot error bars instead of shaded region
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
        options.ErrorBar=false
    end

    % Determine the replicate dimension (the one that is not the time dimension)
    repDim = 3 - dim;  % works for a 2D matrix
    
    % Compute the mean and standard error along the replicate dimension
    m  = mean(data, repDim,"omitmissing");
    se = std(data, 0, repDim,"omitmissing") / sqrt(size(data, repDim));
    
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

    % --- Modified part ---
    if options.ErrorBar
        if isempty(options.DispName_SE)
            se_handle=errorbar(t, m, se, Color=options.Color, LineStyle="none", LineWidth=options.LineWidth);
            se_handle.Annotation.LegendInformation.IconDisplayStyle="off";
        else
            se_handle=errorbar(t, m, se, Color=options.Color, LineStyle="none", LineWidth=options.LineWidth, DisplayName=options.DispName_SE);
        end
    else
        if isempty(options.DispName_SE)
            se_handle=fill([T, fliplr(T)], [Upper, fliplr(Lower)], options.Color, ...
                'EdgeColor', 'none', 'FaceAlpha', options.ShadeAlpha);
            se_handle.Annotation.LegendInformation.IconDisplayStyle="off";
        else
            se_handle=fill([T, fliplr(T)], [Upper, fliplr(Lower)], options.Color, ...
                'EdgeColor', 'none', 'FaceAlpha', options.ShadeAlpha,DisplayName=options.DispName_SE);
        end
    end
    % --- End modification ---

    hold off;

    modes = {'off','on'};
    hold( modes{ holdlog+1 } );

    grid on;
end