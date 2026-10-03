%% Decoding using treebagger, classification, PCA and cPCA
% 
% path2functions='/path/to/matlab_helpers';
% MasterDir="/path/to/data/Test";	% path to store all data
% path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat";
% path2BrainRegions="/path/to/data/BrainRegions.mat";
% p=1;	% set plane for local use
%% TSUBAME
%%%%%%%%%%%%%%%%%%%%%%
path2functions="/path/to/matlab_helpers/";
path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2_For_Tsubame.mat";
path2BrainRegions="/path/to/data/AllMiceData/BrainRegions.mat";
MasterDir="/path/to/results/decode_frontal_parietal/"
% Get Task Index from SGE (or default to 1 if not set)
if isempty(getenv('SGE_TASK_ID'))
    idx = 1;
    disp("Empty SGE_TASK_ID, using idx = 1")
else
    idx = str2double(getenv('SGE_TASK_ID'));
end

p=idx
%%

addpath(genpath(path2functions));
cd(MasterDir);

CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
		dir_BrainRegions_used=path2BrainRegions).PlaneName(p)

%% initialize parfor
% poolobj = gcp('nocreate');
% if isempty(poolobj)
%     poolobj = parpool('Processes');
% end
parpool('local', 12);


%%%%%%%%%%%%%%%%%%%%%
% %% first, obtain the minimum number of cells in each region
% skew_thresh=0.5;
% min_cell_count=Inf;
% for m=1:length(CSP().plane_MouseIND)
% 	temp_log={};
% 	for plane=CSP().plane_MouseIND{m}(:)'
% 		csp=CSP();
% 		csp.load_specified(plane,"Skew","regROI");
% 		skewlog=csp.Skew>=skew_thresh;
% 		alllog={};
% 		for reg=1:length(csp.Reg_inRegLog_reg)
% 			alllog{end+1}=csp.Reg_inRegLog_reg{reg}(:).*skewlog(:);
% 		end
% 		alllog=cellfun(@(a)sum(a,"all"),alllog);
% 		temp_log{end+1}=alllog;
% 	end
% 	temp_log=cat(1,temp_log{:});
% 	temp_log=max(temp_log,[],1);
% 	min_cell_count=min(min_cell_count,min(temp_log));
% end
% min_cell_count

if false
%% first, obtain the minimum number of cells in each region for each plane
skew_thresh=0.5;
cell_count_4_region={};
for m=1:length(CSP().plane_MouseIND)
	for plane=CSP().plane_MouseIND{m}(:)'
		csp=CSP();
		csp.load_specified(plane,"Skew","regROI");
		skewlog=csp.Skew>=skew_thresh;
		alllog={};
		for reg=1:length(csp.Reg_inRegLog_reg)
			alllog{end+1}=csp.Reg_inRegLog_reg{reg}(:).*skewlog(:);
		end
		alllog=cellfun(@(a)sum(a,"all"),alllog);
		cell_count_4_region{end+1}=alllog;
	end
end
cell_count_4_region=cat(1,cell_count_4_region{:});
min_cell_count_4_region=min(cell_count_4_region,[],1)
end
% the minimum # of cells in PPC is 309, and that in S1 trunk is 106

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% actual decoding %%%%%%%%%%%%%%
totaltime=tic;

skew_thresh=0.5
% repetition=5
% numcell2sample=200
groupnum=10
pca_dim2keep=1000
cPCA_weights=[0.3,0.5,0.7,0.9];
cPCA_dim2keep=100
NumTrees=100
NumCPCAaxes2remove=10

currparamDLC=["tDXYZbp_magnitude_mean","eigen_DLC"]
Data2combine="SP only"



cellcount1=[];
cellcount2=[];

i=p;
% for i=1:length(CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
% 		dir_BrainRegions_used=path2BrainRegions).PlaneName)
csp=CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
	dir_BrainRegions_used=path2BrainRegions);
csp.loadData(i, ...
	PeriodCutoff4Derivatives_sec=0.1, ...
	triErrorThresh_mm=4);

csp.setCurrentData("p8");
csp.make_eigen_posture(NumDim2keep=10);
% csp.combine_SP_and_DLC(SP=currparamSP, ...
% 	DLC="");
csp.SP_makeBins(phibincount=2,taubincount=1);
csp.DLC_makeBins(bincount=5,DLC_edge_percent_thresh_For_cutoff=5, ...
	postureFormat="normalPosture",triError_thresh_in_mm=4)
csp.SP_Discretize("theta","phi",EncodeBidirectionalPhi=true);
csp.DLC_Discretize("tXYZbp","tANGbp","tDXYZbp_magnitude","tDXYZbp_magnitude_mean","eigen_DLC");
csp.DLC_Discretize(currparamDLC{:});
csp.combine_Discretized_SP_and_DLC(Data2combine=Data2combine);


%% divide to regions --> PCA --> cPCA
useParallel=true;
region_list={["M1","M2"],["Trunk","Whisker","UL","LL","PPC"]};
rec_data={};
rec_data_base={};

for w=1:numel(cPCA_weights)
	rec_data{w}={};
	for reg=1:numel(region_list)
		ind=get_region_inds(csp,region_list{reg},skew_thresh);
		tckd=csp.tck(:,ind,:,:);
		Csp=copy(csp);
		Csp.split_Data_to_xvalidation_groups(Predictor=tckd,Responser=Csp.Discretized_SP_and_DLC,groupnum=groupnum,dim2split=3);
		Csp.pca_xval_Dataset(Prop2ApplyPCA=["Trainer","Predictor"],Prop2_Multiply_PCA_Axes=["Tester","Predictor"],NumDim2keep=pca_dim2keep);
		Csp.cPCA_xval_Dataset(Prop2ApplycPCA=["Trainer","Predictor_pca"],Prop2_Multiply_cPCA_Axes=["Tester","Predictor_pca"],weight=cPCA_weights(w),NumDim2keep=cPCA_dim2keep);
		%%
		csP=copy(Csp);
		csP.TreeBaggerDecode(NumTrees=NumTrees, ...
			Predictor2Use="Predictor_pca_cpca", ...
			Responser2Use="Responser",...
			UseParallel=useParallel, ...
			Regression=false);
		csP.organize_xval_Dataset(...
			["Tester","Responser"],["Tester","Reconstructed"], ...
			["Trainer","Responser"]);
		rec_data_base{w}{reg}=csP.xval_Dataset_cleared;
		%%
		rec_data{w}{reg}={};
		for Axis=1:NumCPCAaxes2remove
			csP=copy(Csp);
			remove_cPCA_axes(csP,Axis)
			csP.TreeBaggerDecode(NumTrees=NumTrees, ...
				Predictor2Use="Predictor_pca_cpca", ...
				Responser2Use="Responser",...
				UseParallel=useParallel, ...
				Regression=false);
			csP.organize_xval_Dataset(...
				["Tester","Responser"],["Tester","Reconstructed"], ...
				["Trainer","Responser"]);
			rec_data{w}{reg}{Axis}=csP.xval_Dataset_cleared;
		end
		fprintf("%d / %d\n",reg,numel(region_list))
	end
	fprintf("weight %d / %d complete \n",w,numel(cPCA_weights))
end

%% save data

csp.clearData(Properties2keep=["Discretized_SP_and_DLC_names",'Discretized_SP_and_DLC_edges',"currPlaneName4Directory"]);

% save(fullfile(MasterDir,csp.currPlaneName4Directory+".mat"), ...
% 	"rec_data","rec_data_base","cPCA_weights","region_list","csp","-v7.3")
fname = fullfile(MasterDir, csp.currPlaneName4Directory + ".mat");
save(fname, ...        % convert string → char
     'rec_data', ...
     'rec_data_base', ...
     'cPCA_weights', ...
     'region_list', ...
     'csp', ...
     '-v7.3');

toc(totaltime);



%% 
function all_log=get_region_inds(csp,RegNames,skew_thresh)
skew_log=csp.Skew>=skew_thresh;
csP=copy(csp);
Reg_log=csP.Reg_inRegLog_reg(containsSubstring(csP.Reg_RegionNames,RegNames));
if sum(containsSubstring(csP.Reg_RegionNames,RegNames))~=numel(RegNames)
	warning("Region count mismatch")
end
Reg_log=cell2mat(Reg_log(:)');
Reg_log=any(Reg_log,2);
all_log=Reg_log(:)'&skew_log(:)';
end
%%
function remove_cPCA_axes(csp,axis)
for i=1:length(csp.xval_Dataset)
	csp.xval_Dataset{i}.Trainer.Predictor_pca_cpca=remove_axis(csp.xval_Dataset{i}.Trainer.Predictor_pca_cpca,axis);
	csp.xval_Dataset{i}.Tester.Predictor_pca_cpca=remove_axis(csp.xval_Dataset{i}.Tester.Predictor_pca_cpca,axis);
end
end
function newdata=remove_axis(data,axis_num)
IND=makeIND(data,2,setdiff(1:size(data,2),axis_num));
newdata=data(IND{:});
end











