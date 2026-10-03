%% GCARP: movement and posture information in each GCCA component (Fig. 3l, left and right)
% Reads the decoding results of gcarp_decode_movement_posture.m (one file per session).
% For each session and each removed GCCA component (1 to 10), the cross-validated R2 is
% computed relative to the mean of the training responses. The saved value is the R2 with
% all components minus the R2 without the component, for movement (mean speed of the body
% parts, tDXYZbp_magnitude_mean) and for the first eigenposture (eigen_posture_1).
path2functions='/path/to/matlab_helpers';
path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2.mat";
path2BrainRegions="/path/to/data/BrainRegions_with_RSC.mat";
addpath(genpath(path2functions));
OutDir="/path/to/results/gcarp_contributions";

%% load decoding results
folderdir="/path/to/results/gcarp_decode_movement_posture";
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

%% names of the decoded variables (same settings as the decoding job)
currparamSP=["phi"]
currparamDLC=["tDXYZbp_magnitude_mean","eigen_DLC"]

csp=CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
	dir_BrainRegions_used=path2BrainRegions);
csp.loadData(1, ...
	PeriodCutoff4Derivatives_sec=0.1, ...
	triErrorThresh_mm=4);

csp.setCurrentData("p8");
csp.make_eigen_posture(NumDim2keep=10);
csp.combine_SP_and_DLC(SP=currparamSP, ...
	DLC=currparamDLC);

%% assess performance with pseudo Rsquared
param2assess_list=["tDXYZbp_magnitude_mean","eigen_posture_1"];
reglist=["reg_data"];
SubtractFromBaseData=true;

for param2assess=param2assess_list
	for reg=reglist(:)'
		R2_axis_plane=[];
		base_R2_axis_plane=[];
		for plane=1:length(CSP().PlaneName)
			R2_axis=[];
			base_R2_axis=[];
			for Axis=1:numel(alldata{plane}.(reg))
				rec_cat=[];
				resp_cat=[];
				rec_null_cat=[];
				base_rec_cat=[];
				base_resp_cat=[];
				base_rec_null_cat=[];
				for rep=1:numel(alldata{plane}.(reg){Axis})
					for group=1:numel(alldata{plane}.(reg){Axis}{rep})
						ind=csp.SP_and_DLC_combined_names==param2assess;
						Tester=alldata{plane}.(reg){Axis}{rep}{group}.Tester;
						IND=makeIND(Tester.Responser,2,ind);
						resp=alldata{plane}.(reg){Axis}{rep}{group}.Tester.Responser(IND{:});
						rec=alldata{plane}.(reg){Axis}{rep}{group}.Tester.Reconstructed(IND{:});
						if ~isequal(size(resp),size(rec))
							error("Size must match");
						end
						if sum(ind)>1
							error("ind must be one element")
						end

						train_resp=alldata{plane}.(reg){Axis}{rep}{group}.Trainer.Responser(IND{:});
						rec_null=mean(train_resp,[1,3:ndims(train_resp)])+zeros(size(resp));

						% base data
						basereg=reg+"_base";
						baseTester=alldata{plane}.(basereg){rep}{group}.Tester;
						IND=makeIND(baseTester.Responser,2,ind);
						base_resp=baseTester.Responser(IND{:});
						base_rec=baseTester.Reconstructed(IND{:});

						train_base_resp=alldata{plane}.(basereg){rep}{group}.Trainer.Responser(IND{:});
						base_rec_null=mean(train_base_resp,[1,3:ndims(train_base_resp)])+zeros(size(base_resp));

						resp_cat=cat(1,resp_cat,resp(:));
						rec_cat=cat(1,rec_cat,rec(:));
						rec_null_cat=cat(1,rec_null_cat,rec_null(:));
						base_resp_cat=cat(1,base_resp_cat,base_resp(:));
						base_rec_cat=cat(1,base_rec_cat,base_rec(:));
						base_rec_null_cat=cat(1,base_rec_null_cat,base_rec_null(:));
					end
				end
				R2_axis(Axis)=1-(sum((rec_cat-resp_cat).^2)./sum((rec_null_cat-resp_cat).^2));
				base_R2_axis(Axis)=1-(sum((base_rec_cat-base_resp_cat).^2)./sum((base_rec_null_cat-base_resp_cat).^2));
			end
			R2_axis_plane(:,plane)=R2_axis(:);
			base_R2_axis_plane(:,plane)=base_R2_axis(:);
		end
		if SubtractFromBaseData
			R2_axis_plane=base_R2_axis_plane-R2_axis_plane;
		end
	end
	% save data
	delta_R2_axis_plane=R2_axis_plane;
	save(fullfile(OutDir,param2assess+".mat"),"delta_R2_axis_plane")
end
