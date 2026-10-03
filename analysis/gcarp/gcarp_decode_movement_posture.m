%% Decoding using treebagger, regression, anchor of GCCA
% 
% path2functions='/path/to/matlab_helpers';
% MasterDir="/path/to/data/Test";	% path to store all data
% path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat";
% path2BrainRegions="/path/to/data/BrainRegions_with_RSC.mat";
% p=1;	% set plane for local use
%% TSUBAME
%%%%%%%%%%%%%%%%%%%%%%
path2functions="/path/to/matlab_helpers/";
path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2_For_Tsubame.mat";
path2BrainRegions="/path/to/data/AllMiceData/BrainRegions_with_RSC.mat";
MasterDir="/path/to/results/gcarp_decode_movement_posture/"
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
% sample 500 cells if possible, pca to 100 dimensions, and cca 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% actual decoding %%%%%%%%%%%%%%
totaltime=tic;

skew_thresh=0.5
repetition=5
numcell2sample=200
groupnum=10
pca_dim2keep=100
NumTrees=100
NumCCAaxes2remove=10


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
%%
currparamSP=["phi"]
currparamDLC=["tDXYZbp_magnitude_mean","eigen_DLC"]
currsetting="SP_and_tDmag_eig"
params2decode=["tXYZbp","tANGbp","tDXYZbp_magnitude","tDXYZbp_magnitude_mean","eigen_DLC"]
% Data2combine="Both"
% Data2combinelist=["Both","SP only","DLC only"];

% currparam={"tXYZbp"}
% currsetting="tXYZbp"
% params2decode=["tXYZbp","tANGbp","tDXYZbp_magnitude","tDXYZbp_magnitude_mean","eigen_DLC"];
% Data2combine="DLC only"
% Data2combinelist=["Both","SP only","DLC only"];


csp.setCurrentData("p8");
csp.make_eigen_posture(NumDim2keep=10);
csp.combine_SP_and_DLC(SP=currparamSP, ...
	DLC=currparamDLC);
% csp.SP_makeBins(phibincount=5,taubincount=1);
% csp.DLC_makeBins(bincount=5,DLC_edge_percent_thresh_For_cutoff=5, ...
% 	postureFormat="normalPosture",triError_thresh_in_mm=4)
% csp.SP_Discretize("theta","phi",EncodeBidirectionalPhi=false);
% % csp.DLC_Discretize("tXYZbp","tANGbp","tDXYZbp_magnitude","tDXYZbp_magnitude_mean","eigen_DLC");
% csp.DLC_Discretize(currparam{:});
% csp.combine_Discretized_SP_and_DLC(Data2combine=Data2combine);

%% repeat for axis and rep
useParallel=true;
S1trunk_data_base={};
PPC_data_base={};
S1trunk_data={};
PPC_data={};
cellcount=[];

listOfRegions={"PPC","M1","M2","Trunk","Whisker","UL","LL","RSC"};

regnames=cellfun(@(a)join(a,"_"),listOfRegions,UniformOutput=true);

% initialize arrays
reg_data={};
reg_data_base={};

for Axis=0:NumCCAaxes2remove
	axtime=tic;
	csp_reg_reparray=CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
		dir_BrainRegions_used=path2BrainRegions).empty(0,repetition);
	for rep=1:repetition
		curr_regnames=regnames;

		reptime=tic;
		csp_reg=struct;
		for s_ind=1:numel(listOfRegions)
			s=curr_regnames(s_ind);
			[csp_reg.(s),~]=make_DataSets_forRegions(csp,listOfRegions{s_ind}, ...
				skew_thresh,rep,numcell2sample,groupnum,pca_dim2keep);
		end
		valid_reg_log=[];
		for s=curr_regnames(:)'
			valid_reg_log(end+1)=~isempty(csp_reg.(s));
			if isempty(csp_reg.(s))
				csp_reg=rmfield(csp_reg,s);		% only the regions with enough # of cells
			end
		end
		curr_regnames=regnames(logical(valid_reg_log));	% only the regions with enough # of cells

		csp.gcca_xval_Dataset(csp_reg,...
			Prop2ApplyCCA=["Trainer","Predictor_pca"], ...
			Prop2_Multiply_CCA_Axes=["Tester","Predictor_pca"], ...
			ObtainAnchor=true);
		structfun(@(csp)remove_CCA_axis(csp,Axis),csp_reg,UniformOutput=false);
		for s=curr_regnames(1)
			fprintf("%s\n",s);
			csP=csp_reg.(s);
			csP.TreeBaggerDecode(NumTrees=NumTrees, ...
				Predictor2Use="Predictor_pca_ccaAnchor", ...
				Responser2Use="Responser",...
				UseParallel=useParallel, ...
				Regression=true);
			csP.organize_xval_Dataset(...
				["Tester","Responser"],["Tester","Reconstructed"], ...
				["Tester","Predictor_pca_ccaCorr"],["Trainer","Predictor_pca_ccaCorr"], ...
				["Trainer","Responser"]);
			csp_reg_reparray(rep)=copy(csP);
		end
		disp("%%%%%%%%%%%%%")
		fprintf("Rep %d / %d complete \n", ...
			rep,repetition);
		toc(reptime)
		disp("%%%%%%%%%%%%%")
	end
	if Axis==0
		reg_data_base={csp_reg_reparray.xval_Dataset_cleared};
	else
		reg_data{Axis}={csp_reg_reparray.xval_Dataset_cleared};
	end

	fprintf("Axes %d / %d complete ################### \n", ...
			Axis,NumCCAaxes2remove);
	toc(axtime);
	disp("########################################")
end
% end

%% save data
% save decoded data
% if Data2combine=="SP only"
% 	currsetting="SP";
% else
% 	currsetting=currparam;
% end

save(fullfile(MasterDir,csp.currPlaneName4Directory+"_"+currsetting+".mat"), ...
	"reg_data_base","reg_data","curr_regnames","-v7.3")

% save additional data (csp)
csp.clearAndSave(...
	Properties2keep=["Discretized_SP_and_DLC_names", ...
	"Discretized_SP_and_DLC_undiscretized_names",...
	"Discretized_SP_and_DLC_edges"], ...
	filename=fullfile(MasterDir,csp.currPlaneName4Directory+"_"+currsetting+"_"+"info.mat"), ...
	SaveAs_v7p3=true)

toc(totaltime);


%% function to remove axis from cca
function remove_CCA_axis(csp,axis)
arguments
	csp
	axis=0
end
if axis==0
	return;
end
for i=1:numel(csp.xval_Dataset)
	data=csp.xval_Dataset{i}.Trainer.Predictor_pca_cca;
	IND=makeIND(data,2,setdiff(1:size(data,2),axis));
	newdata=data(IND{:});
	csp.xval_Dataset{i}.Trainer.Predictor_pca_cca=newdata;

	data=csp.xval_Dataset{i}.Trainer.Predictor_pca_ccaAxes;
	IND=makeIND(data,2,setdiff(1:size(data,2),axis));
	newdata=data(IND{:});
	csp.xval_Dataset{i}.Trainer.Predictor_pca_ccaAxes=newdata;

	data=csp.xval_Dataset{i}.Trainer.Predictor_pca_ccaAnchor;
	IND=makeIND(data,2,setdiff(1:size(data,2),axis));
	newdata=data(IND{:});
	csp.xval_Dataset{i}.Trainer.Predictor_pca_ccaAnchor=newdata;

	data=csp.xval_Dataset{i}.Tester.Predictor_pca_cca;
	IND=makeIND(data,2,setdiff(1:size(data,2),axis));
	newdata=data(IND{:});
	csp.xval_Dataset{i}.Tester.Predictor_pca_cca=newdata;

	data=csp.xval_Dataset{i}.Tester.Predictor_pca_ccaAnchor;
	IND=makeIND(data,2,setdiff(1:size(data,2),axis));
	newdata=data(IND{:});
	csp.xval_Dataset{i}.Tester.Predictor_pca_ccaAnchor=newdata;
end
end

%% function to split data
function [csP,numcell]=make_DataSets_forRegions(csp,RegNames, ...
	skew_thresh,rep,numcell2sample,groupnum,pca_dim2keep)
% initialize rng
rng(rep);

skew_log=csp.Skew>=skew_thresh;
csP=copy(csp);
Reg_log=csP.Reg_inRegLog_reg(containsSubstring(csP.Reg_RegionNames,RegNames));
if sum(containsSubstring(csP.Reg_RegionNames,RegNames))~=numel(RegNames)
	warning("Region count mismatch")
end
Reg_log=cell2mat(Reg_log(:)');
Reg_log=any(Reg_log,2);
all_log=Reg_log(:)'&skew_log(:)';
all_ind=find(all_log);
if numel(all_ind)<numcell2sample
	csP=[];
	numcell=[];
	return;
end
numcell2sample=min(numcell2sample,numel(all_ind));	% the max # of cells to sample are limited to the # of cells inside the region
sampleIND=randperm(numel(all_ind),numcell2sample);
all_ind=all_ind(sampleIND);
IND=makeIND(csP.tck,2,all_ind);
csP.tck=csP.tck(IND{:});

csP.split_Data_to_xvalidation_groups(...
	Predictor=csP.tck, ...
	Responser=csp.SP_and_DLC_combined, ...
	groupnum=groupnum);

csP.pca_xval_Dataset( ...
	Prop2ApplyPCA=["Trainer","Predictor"], ...
	Prop2_Multiply_PCA_Axes=["Tester","Predictor"],...
	NumDim2keep=pca_dim2keep);

numcell=sum(all_log);
end












