%% Frontal and parietal tilt-direction decoding (Fig. 3f-i, Extended Data Fig. 10k)
% Reads the decoding results of TB_cPCA_20250706 (one file per session). For each session
% and cPCA weight, collects the decoded and true tilt directions (radians) of the decoders
% that use all 100 cPCA dimensions, for M1 and M2 (frontal) and for S1 and PPC (parietal),
% at time bins 10 to 20 of every test trial. The figure workflows read Reg_comparison.mat.
path2functions='/path/to/matlab_helpers';
addpath(genpath(path2functions));
OutDir="/path/to/data/TB_classi_cPCA_20250706";

%% load decoding results
folderdir="/path/to/data/TB_classi_cPCA_20250707";
filelist=getAllFiles(folderdir);

alldata={};

for i=1:length(CSP().PlaneName)
	plane=CSP().PlaneName4Directory(i);
	ind=containsSubstring(filelist,plane);
	if sum(ind)>2
		error("too many files")
	end
	temp=load(filelist{ind});
	alldata{i}=temp;
	i
end

csp=copy(temp.csp);

%% decoded and true directions of the frontal and parietal decoders
intra_trial_ind=10:20;
regnames=["M1_M2", "PPC_S1"];
param2assess="theta";	% corresponds to 1 parameter

recdata={};
respdata={};
for reg2plot=1:2
	recdata.(regnames(reg2plot))={};
	respdata.(regnames(reg2plot))={};
	for plane=1:length(CSP().PlaneName)
		for wIND_cPCA=1:numel(alldata{plane}.cPCA_weights)
			base_rec_cat=[];
			base_resp_cat=[];
			for g=1:numel(alldata{plane}.rec_data_base{wIND_cPCA}{reg2plot})
				ind=contains(csp.Discretized_SP_and_DLC_names,param2assess);
				edges=csp.Discretized_SP_and_DLC_edges{ind}{1};
				f=@(a)Undiscretize(round(a),edges);
				base_rec=alldata{plane}.rec_data_base{wIND_cPCA}{reg2plot}{g}.Tester.Reconstructed(intra_trial_ind,ind,:,:);
				base_resp=alldata{plane}.rec_data_base{wIND_cPCA}{reg2plot}{g}.Tester.Responser(intra_trial_ind,ind,:,:);
				base_rec=f(base_rec);base_resp=f(base_resp);

				base_rec_cat=cat(1,base_rec_cat,base_rec(:));
				base_resp_cat=cat(1,base_resp_cat,base_resp(:));
			end
			recdata.(regnames(reg2plot)){plane,wIND_cPCA}=base_rec_cat;
			respdata.(regnames(reg2plot)){plane,wIND_cPCA}=base_resp_cat;
		end
	end
end

%% save data
cPCA_weights=alldata{plane}.cPCA_weights;
save(fullfile(OutDir,"Reg_comparison.mat"),"recdata","respdata","cPCA_weights")
