function mask=generate_mask4conv(x_pix_per_mm,y_pix_per_mm,radius_in_um,option,padding_in_mm)
% the resulting mask has an volume of 1 mm^3
% --> convolute with 2d histogram of cell counts per pixel^2 and get cell
% counts per mm^2
arguments
	x_pix_per_mm
	y_pix_per_mm
	radius_in_um
	option (1,1) string {mustBeMember(option, ["Gaussian","Ellipse"])}
	padding_in_mm=0
end
radius_in_mm=radius_in_um/1000;

radx=round(radius_in_mm*x_pix_per_mm);
rady=round(radius_in_mm*y_pix_per_mm);

pad=padding_in_mm;
Radx=round(radx+pad*x_pix_per_mm);
Rady=round(rady+pad*y_pix_per_mm);

elps=[1:2*Radx-1]+1i.*[1:2*Rady-1]';
switch option
	case "Ellipse"
		mask=((1/radx^2)*real(elps-(Radx+1i*Rady)).^2+(1/rady^2)*imag(elps-(Radx+1i*Rady)).^2)<=1;
		% mask=mask./(pi*radx*rady);	% 1 pix^2
	case "Gaussian"
		cent=Radx+Rady*1i;
		mask=(1/(2*pi*radx*rady))*exp(-0.5*((real(elps-cent)/radx).^2+(imag(elps-cent)/rady).^2));	% integrate by pix--> 1 --> cell counts / pix^2
end
% mask=mask*x_pix_per_mm*y_pix_per_mm;	% convert from 1pix^2 to 1mm^2



end

%%
function testThisFunc
%%
figure;
imagesc(mask)
end