function newtpk=PreprocessSPdata_Accurate(org)
tpk=org.sp_t;
% 2nd dimension of tpk must be the parameters
IND=repmat({':'},1,ndims(tpk));

IND{2}=1;
theta=tpk(IND{:});
IND{2}=2;
phi=tpk(IND{:});
theta=atan2(org.sinsp_t(:,1),org.cossp_t(:,1));


theta(phi<0)=theta(phi<0)+pi;
phi(phi<0)=-phi(phi<0);

% Process theta (limit the range)
theta=ThetaMod4DiscretizationUsingmakeThetaBins(theta);


IND{2}=1;
tpk(IND{:})=theta;
IND{2}=2;
tpk(IND{:})=phi;
newtpk=tpk;


end