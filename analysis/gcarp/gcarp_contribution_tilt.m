%% GCARP: ground-tilt information in each GCCA component (Fig. 3l, middle)
% Reads the decoding results of gcarp_decode_tilt.m (one file per session).
% For each session and each removed GCCA component (1 to 10), the score is the fraction
% of test samples whose decoded tilt direction is within 45 deg of the true direction.
% The saved value is the score with all components minus the score without the component.
path2functions='/path/to/matlab_helpers';
addpath(genpath(path2functions));
OutDir="/path/to/results/gcarp_contributions";

%% load decoding results
folderdir="/path/to/results/gcarp_decode_tilt";
filelist=getAllFiles(folderdir);

alldata={};
allinfo={};

for i=1:length(CSP().PlaneName)
	plane=CSP().PlaneName4Directory(i);
	ind=containsSubstring(filelist,plane+"_SP_and_tDmag_eig");
	altind=containsSubstring(filelist,plane+"_SP_and_tDmag_eig_info");
	if sum(ind)>2
		error("too many files")
	end
	dataind=logical(ind-altind);
	if sum(dataind)>1
		error("too many files")
	end
	temp=load(filelist{dataind});
	alldata{i}=temp;
	temp=load(filelist{altind});
	allinfo{i}=temp.obj;
	i
end

%% assess performance with delta deg, theta
deg_thresh=45;
ddeg_thresh=10;	% tolerance (keep this low)
trial_time2use=10:15;

rad_thresh=deg2rad(deg_thresh);
drad_thresh=deg2rad(ddeg_thresh);

param2assess="theta";	% corresponds to 1 parameter
reglist=["reg_data"];
SubtractFromBaseData=true;

for reg=reglist(:)'
	ddegScore_axis_plane=[];
	base_ddegScore_axis_plane=[];
	for plane=1:length(CSP().PlaneName)
		ddegScore_axis=[];
		base_ddegScore_axis=[];
		csp=allinfo{plane};
		for Axis=1:numel(alldata{plane}.(reg))
			rec_cat=[];
			resp_cat=[];
			base_rec_cat=[];
			base_resp_cat=[];
			for rep=1:numel(alldata{plane}.(reg){Axis})
				for group=1:numel(alldata{plane}.(reg){Axis}{rep})
					ind=csp.Discretized_SP_and_DLC_names==param2assess;
					Tester=alldata{plane}.(reg){Axis}{rep}{group}.Tester;
					IND=makeIND(Tester.Responser,2,ind,1,trial_time2use);
					resp=alldata{plane}.(reg){Axis}{rep}{group}.Tester.Responser(IND{:});
					rec=alldata{plane}.(reg){Axis}{rep}{group}.Tester.Reconstructed(IND{:});
					if ~isequal(size(resp),size(rec))
						error("Size must match");
					end
					if sum(ind)>1
						error("ind must be one element")
					end

					% base data
					basereg=reg+"_base";
					baseTester=alldata{plane}.(basereg){rep}{group}.Tester;
					IND=makeIND(baseTester.Responser,2,ind,1,trial_time2use);
					base_resp=baseTester.Responser(IND{:});
					base_rec=baseTester.Reconstructed(IND{:});

					resp_cat=cat(1,resp_cat,resp(:));
					rec_cat=cat(1,rec_cat,rec(:));
					base_resp_cat=cat(1,base_resp_cat,base_resp(:));
					base_rec_cat=cat(1,base_rec_cat,base_rec(:));
				end
			end
			% undiscretize variables
			uresp_cat=Undiscretize(resp_cat,allinfo{plane}.Discretized_SP_and_DLC_edges{ind}{1});
			urec_cat=Undiscretize(rec_cat,allinfo{plane}.Discretized_SP_and_DLC_edges{ind}{1});
			base_uresp_cat=Undiscretize(base_resp_cat,allinfo{plane}.Discretized_SP_and_DLC_edges{ind}{1});
			base_urec_cat=Undiscretize(base_rec_cat,allinfo{plane}.Discretized_SP_and_DLC_edges{ind}{1});

			ddeg=wrapToPi(urec_cat-uresp_cat);
			ddegScore_axis(Axis)=sum(abs(ddeg)<=(rad_thresh+drad_thresh))/numel(ddeg);

			ddeg=wrapToPi(base_urec_cat-base_uresp_cat);
			base_ddegScore_axis(Axis)=sum(abs(ddeg)<=(rad_thresh+drad_thresh))/numel(ddeg);
		end
		ddegScore_axis_plane(:,plane)=ddegScore_axis(:);
		base_ddegScore_axis_plane(:,plane)=base_ddegScore_axis(:);
	end
	if SubtractFromBaseData
		ddegScore_axis_plane=base_ddegScore_axis_plane-ddegScore_axis_plane;
	end
end
%% save data
delta_degScore_axis_plane=ddegScore_axis_plane;
save(fullfile(OutDir,"theta.mat"),"delta_degScore_axis_plane");
