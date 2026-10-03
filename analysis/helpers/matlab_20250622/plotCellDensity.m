function plotCellDensity(regROI, cellLog, mask4conv, options)
% plotCellDensity overlays a smoothed density ratio map of "qualified" cells 
% on top of a background image (if provided) or a blank canvas.
%
% This function performs the following steps:
%   1. Bins the cell coordinates (given in regROI) into pixels using histcounts2.
%   2. Computes two binned maps: one for all cells and one for qualified cells 
%      (as indicated by the logical vector cellLog).
%   3. Smooths these binned maps by convolving them with a provided mask (mask4conv),
%      for example, a 2D Gaussian kernel.
%   4. Computes a ratio map (qualified-cell density / overall cell density) for each pixel.
%   5. Optionally applies a threshold to the overall density before displaying the ratio.
%   6. Displays the density ratio overlay using imagesc (with a user-specified colormap
%      and transparency) on top of the background image (if provided).
%
% USAGE:
%   plotCellDensity(regROI, cellLog, mask4conv, options)
%
% Required Inputs:
%   regROI   - An N-by-2 numeric array of cell coordinates. Column 1 contains the 
%              y-coordinates (vertical; top at 1) and column 2 contains the x-coordinates.
%   cellLog  - An N-by-1 logical vector indicating the "qualified" cells.
%   mask4conv- A 2D numeric array used as a convolution kernel (e.g., a Gaussian kernel)
%              to smooth the binned cell count maps.
%
% Name–Value Pair Options (provided via the options structure):
%   options.AtlasImage      - (Default: []) Numeric image array to be used as a background.
%                             If empty, a blank background is used.
%
%   options.overlayColormap - (Default: jet(256)) Numeric colormap array for the overlay.
%
%   options.alphaVal        - (Default: 1) Scalar in [0,1] specifying the transparency
%                             of the overlay (1 = fully opaque).
%
%   options.mask_thresh     - (Default: []) Numeric threshold. Only pixels in the overall 
%                             density map with values greater than or equal to this value
%                             are retained in the ratio map.
%
%   options.Title           - (Default: '') Character vector to set the plot title.
%
%   options.imagescmax      - (Default: 1) Numeric value setting the maximum value for
%                             scaling the imagesc overlay.
%

% Define input arguments and their default values using an arguments block.
arguments
    regROI (:,2) double
    cellLog (:,1) logical
    mask4conv (:,:) double
    options.AtlasImage {mustBeNumeric} = []      
	options.AtlasImage_mask=[]
    options.overlayColormap {mustBeNumeric} = jet(256)       
    options.alphaVal (1,1) double {mustBeGreaterThanOrEqual(options.alphaVal,0), ...
                          mustBeLessThanOrEqual(options.alphaVal,1)} = 1  
    options.mask_thresh {mustBeNumeric} = []               
    options.Title (1,:) char = ''                          
    options.imagescmax (1,1) double = 1                    
	options.Plot_mask4conv (1,1) logical=true
end

% Assign the options to internal variables for clarity.
AtlasImage      = options.AtlasImage;
overlayColormap = options.overlayColormap;
alphaVal        = options.alphaVal;
mask_thresh     = options.mask_thresh;
Title           = options.Title;
imagescmax      = options.imagescmax;

%% Set up the binning based on the background image (if provided)
if ~isempty(AtlasImage)
    % Get the image size for binning.
    [nRows, nCols, ~] = size(AtlasImage);
    % Define bin edges so that each pixel (with center at an integer value) gets its own bin.
    yEdges = 0.5:1:(nRows+0.5);
    xEdges = 0.5:1:(nCols+0.5);
    
    % Use the supplied alpha value for the overlay.
    alphaUsed = alphaVal;
    
    combim = AtlasImage;  %#ok<NASGU> (for clarity, this variable holds the background image)
else
    % No AtlasImage provided; determine bin edges from the maximum cell coordinates.
    maxY = ceil(max(regROI(:,1)));
    maxX = ceil(max(regROI(:,2)));
    yEdges = 0.5:1:(maxY+0.5);
    xEdges = 0.5:1:(maxX+0.5);
    
    % Use full opacity when there is no background image.
    alphaUsed = 1;
    
    combim = zeros(maxY, maxX);  %#ok<NASGU>
end

%% Compute binned counts and smoothed densities

% Binning for all cells:
[allCounts,~,~] = histcounts2(regROI(:,1), regROI(:,2), yEdges, xEdges);
roidense = conv2(allCounts, mask4conv, 'same');  % Smoothed density of all cells

% Binning for qualified cells:
[qualCounts,~,~] = histcounts2(regROI(cellLog,1), regROI(cellLog,2), yEdges, xEdges);
roipercent = conv2(qualCounts, mask4conv, 'same');  % Smoothed density of qualified cells

% Compute the ratio (qualified density / all-cell density)
roipercent = roipercent ./ roidense;
roipercent(isnan(roipercent)) = 0;  % Replace NaNs (from 0/0) with zeros

% Apply mask threshold if provided.
if ~isempty(mask_thresh)
    mask = roidense >= mask_thresh;
    roipercent = roipercent .* double(mask);
end

%% Display the background image and overlay the density ratio map


if ~options.Plot_mask4conv
	if ~isempty(AtlasImage)
		% If a background image is provided, display it.
		atlas=imshow(AtlasImage);
		hold on;
		if ~isempty(options.AtlasImage_mask)
			mask=zeros(size(AtlasImage,[1,2])+[0,size(mask4conv,2)]);
			mask(:,size(mask4conv,2)+1:end)=options.AtlasImage_mask;
			atlas.AlphaData=mask*alphaUsed;
		end
	end

	% Plot the density ratio overlay using imagesc.
	h = imagesc(roipercent, [0, imagescmax]);
	axis equal;
	if isempty(options.AtlasImage_mask)
		h.AlphaData = alphaUsed;  % Set overlay transparency.
	end
	title(Title);
	colormap(gca(),overlayColormap);
	colorbar;
	hold off;
else
	if ~isempty(AtlasImage) && isempty(options.AtlasImage_mask)
		% If a background image is provided, display it.
		I=zeros(size(AtlasImage)+[0,size(mask4conv,2),0],"uint8");
		I(:,size(mask4conv,2)+1:end,:)=AtlasImage;
		atlas=imshow(I);
		hold on;
		if ~isempty(options.AtlasImage_mask)
			mask=zeros(size(AtlasImage,[1,2])+[0,size(mask4conv,2)]);
			mask(:,size(mask4conv,2)+1:end)=options.AtlasImage_mask;
			atlas.AlphaData=mask*alphaUsed;
		end
	end
	im=zeros([size(roipercent,1),size(roipercent,2)+size(mask4conv,2)]);
	im(:,size(mask4conv,2)+1:end)=roipercent;
	% im(1:size(mask4conv,1),1:size(mask4conv,2))=mask4conv;

	% Plot the density ratio overlay using imagesc.
	h = imagesc(im, [0, imagescmax]);
	if isempty(options.AtlasImage_mask)
		h.AlphaData = alphaUsed;  % Set overlay transparency.
	elseif ~isempty(AtlasImage)
		% show atlas image overlay with transparency
		I=zeros(size(AtlasImage)+[0,size(mask4conv,2),0],"uint8");
		I(:,size(mask4conv,2)+1:end,:)=AtlasImage;
		hold on;
		atlas=imshow(I);
		if ~isempty(options.AtlasImage_mask)
			mask=zeros(size(AtlasImage,[1,2])+[0,size(mask4conv,2)]);
			mask(:,size(mask4conv,2)+1:end)=options.AtlasImage_mask;
			atlas.AlphaData=mask*alphaUsed;
		end
	end
	axis equal;
	title(Title);
	colormap(gca(),overlayColormap);
	colorbar;
	imagesc(mask4conv)
	hold off;
end
impixelinfo;

end
