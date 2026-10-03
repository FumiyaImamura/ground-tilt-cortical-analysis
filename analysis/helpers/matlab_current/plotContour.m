function Fig=plotContour(options)
arguments
	options.contour
	options.img (:,:,:) {mustBeNumeric}=[]	% usually AtlasImage
	options.roi (:,:) {mustBeNumeric}=[]
end
contour=options.contour;

Fig=figure;
if ~isempty(options.img)
	img=options.img;
	if size(img,3)>1
		imshow(img)
	else
		lowprcnt=prctile(img,15,'all');
		highprcnt=prctile(img,99,'all');
		img(img<lowprcnt)=lowprcnt;
		img(img>highprcnt)=highprcnt;
		img=double(img-lowprcnt)./double(highprcnt-lowprcnt).*double(intmax("uint16"));
		img=uint16(img);
	end
	hold on;
end

if iscell(contour)
	for i=1:length(contour)
		[B, L] = bwboundaries(contour{i}, 'noholes');
		boundary=B{1};
		plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
		hold on;
	end
else
	[B, L] = bwboundaries(contour, 'noholes');
	boundary=B{1};
	plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2);
	hold on;
end

if ~isempty(options.roi)
	roi=options.roi;
	plot(roi(:,2), roi(:,1), 'r.');
end


end