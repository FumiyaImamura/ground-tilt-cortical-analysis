function inRegion_g=assessInsideRegion_g(ROI,Contour_g)
% inRegion: 1 x g: logical vector representing if a cell is inside the
% region or not.

arguments
	ROI (:,:) {mustBeNumeric}		% same as indexing in Atlas_2p.Coord_IND
	Contour_g (:,1) {mustBeA(Contour_g,"cell")}
end

% reshape ROI so that n x 2
if size(ROI,1)==2
	ROI=ROI';
elseif size(ROI,1)~=2 && size(ROI,2)~=2
	error('ROI must contain a dimension with a size 2');
end

if ~isvector(Contour_g)
	error('Contour_g must be a cell vector')
end

inRegion_g=cell(length(Contour_g),1);
for i=1:length(Contour_g)
	a=assessInsideRegion(ROI,Contour_g{i});
	inRegion_g{i}=a(:);
end

end