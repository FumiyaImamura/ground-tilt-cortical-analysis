function updated_theta_array=ThetaMod4DiscretizationUsingmakeThetaBins(original_theta_array)
% original_theta_array refers to an array of any dimension corresponding to
% theta
% For example, tpk(:,1,:)
% since I want it to span from pi/8 to pi/8+2*pi, this function calculates
% the mod of the original theta array to match this range

% create thetabins
thetaedges=makeThetaBins();
a=thetaedges(1);
updated_theta_array=mod(original_theta_array-a,2*pi)+a;

end