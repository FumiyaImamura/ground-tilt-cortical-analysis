function inRegion=assessInsideRegion(ROI,Contour)
% inRegion: n x 1: logical vector representing if a cell is inside the
% region or not

arguments
	ROI (:,:) {mustBeNumeric}		% same as indexing in Atlas_2p.Coord_IND
	Contour (:,:) {mustBeA(Contour,"logical")}
end

% reshape ROI so that n x 2
if size(ROI,1)==2
	ROI=ROI';
elseif size(ROI,1)~=2 && size(ROI,2)~=2
	error('ROI must contain a dimension with a size 2');
end


[B, L] = bwboundaries(Contour, 'noholes');

boundary=B{1};
if length(B)>1
	error('There are multiple boundaries')
end

inRegion = inpolygon(ROI(:,2),ROI(:,1), boundary(:,2), boundary(:,1));

end