%% Ridge regression of continuous data

% path2functions='/path/to/matlab_helpers';
% MasterDir="/path/to/data/Test";	% path to store all data
% path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2.mat";
% path2BrainRegions="/path/to/data/BrainRegions.mat";
% p=1;	% set plane for local use
%% TSUBAME
%%%%%%%%%%%%%%%%%%%%%%
path2functions="/path/to/matlab_helpers/";
path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2_For_Tsubame.mat";
path2BrainRegions="/path/to/data/AllMiceData/BrainRegions.mat";
MasterDir="/path/to/results/ridge_eigenposture_ablation/"
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
% parpool('local', 12);


%% ridge regression
ticmaster=tic;


sigma=1;
% ridgeparams=10.^[-12:2:12];
ridgeparams=10.^[2:0.5:8];
groupnum=10;
make_Data_single=false;
numcellac2saveAsAnExample=10;

csp=CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
	dir_BrainRegions_used=path2BrainRegions);

% file2saveData=fullfile(MasterDir,append(csp.PlaneName4Directory(p),".mat"));

%%
i=p;

% mkdir(folder2saveData);
tstart=tic;
% for p=1:length(csp.PlaneName)
tic
csp.loadData(i, ...
	PeriodCutoff4Derivatives_sec=0.1, ...
	triErrorThresh_mm=4);




%% ridge using all parameters
fprintf("starting ridge using all parameters...\n")
tstart=tic;
csp.loadData_continuous_time(UseCachedData=true);
csp.setCurrentData_to_continuous()
csp.SP_makeBins(phibincount=5,taubincount=5,Zbincount=5)
csp.SP_Discretize(["theta","phi"],"tau","Z",EncodeBidirectionalPhi=true,EncodeBidirectional=true);
csp.SP_GaussianHotEncode_Discretized(sigma)
csp.DLC_makeBins(bincount=5,DLC_edge_percent_thresh_For_cutoff=5, ...
	postureFormat="normalPosture",triError_thresh_in_mm=4)
csp.DLC_Discretize("tXYZbp","tDXYZbp","tD2XYZbp","tANGbp","tDANGbp","tD2ANGbp", ...
	"tDXYZbp_magnitude","tDANGbp_magnitude","tD2XYZbp_magnitude","tD2ANGbp_magnitude")
csp.DLC_GaussianHotEncode_Discretized(sigma)
csp.combinePredictors("Both",true,removeSmilarParameters=false)

csp.split_Data_to_xvalidation_groups_continuous( ...
	Predictor=csp.all_GaussianHot, ...
	Responser=csp.cellac_continuous, ... 
	groupnum=groupnum, ...
	TimeChunkSec=5);
csp.RidgeDecode(...
	Predictor2Use="Predictor", ...
	Responser2Use="Responser", ...
	make_Data_single=make_Data_single, ...
	ridgeparams=ridgeparams, ...
	UseGPU=false);
Both_score_ridge_cells=csp.xval_assess_performance_ridge(ReferenceResponser="Responser");

% save cellac with high scores
csp2save=saveExampleCellac(csp,Both_score_ridge_cells.R,numcellac2saveAsAnExample);

toc(tstart);


%%
csp=CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
	dir_BrainRegions_used=path2BrainRegions);
csp.loadData(p, ...
	PeriodCutoff4Derivatives_sec=0.1, ...
	triErrorThresh_mm=4, ...
	saveRawQ_bp=true);
csP=copy(csp);
%% ridge using dropout eigenpostures
eigenpostureDim2analyze=10;
score_eigdim=[];
for dim=1:eigenpostureDim2analyze
	fprintf("starting ridge using eigenpostures with no dim %d...\n",dim)
	tstart=tic;
	csp=copy(csP);
	csp.loadData_continuous_time(UseCachedData=true);
	csp.setCurrentData_to_continuous()
	Csp=copy(csp);
	Csp.make_eigen_posture();
	posture=Csp.reconstruct_from_eigen_posture(PostureData=csp.RawQ_bp,Dims2use=setdiff(1:42,dim));
	DLC_raw=csp.make_modified_DLC_raw(posture,csp.rec3d);
	csp.loadData_continuous_time(UseCachedData=true,Cached_DLC_raw=DLC_raw);
	csp.setCurrentData_to_continuous()

	csp.SP_makeBins(phibincount=5,taubincount=5,Zbincount=5)
	csp.SP_Discretize(["theta","phi"],"tau","Z",EncodeBidirectionalPhi=true,EncodeBidirectional=true);
	csp.SP_GaussianHotEncode_Discretized(sigma)
	csp.DLC_makeBins(bincount=5,DLC_edge_percent_thresh_For_cutoff=5, ...
		postureFormat="normalPosture",triError_thresh_in_mm=4)
	csp.DLC_Discretize("tXYZbp","tDXYZbp","tD2XYZbp","tANGbp","tDANGbp","tD2ANGbp", ...
		"tDXYZbp_magnitude","tDANGbp_magnitude","tD2XYZbp_magnitude","tD2ANGbp_magnitude")
	csp.DLC_GaussianHotEncode_Discretized(sigma)
	csp.combinePredictors("Both",true,removeSmilarParameters=false)

	csp.split_Data_to_xvalidation_groups_continuous( ...
		Predictor=csp.all_GaussianHot, ...
		Responser=csp.cellac_continuous, ...
		groupnum=groupnum, ...
		TimeChunkSec=5);
	csp.RidgeDecode(...
		Predictor2Use="Predictor", ...
		Responser2Use="Responser", ...
		make_Data_single=make_Data_single, ...
		ridgeparams=ridgeparams, ...
		UseGPU=false);
	score_eigdim{dim}=csp.xval_assess_performance_ridge(ReferenceResponser="Responser");

	% figure;
	% scatax=[];
	% plotax=[];
	% for i=1:100:size(posture,1)
	% 	scatax=csp.scatter3_posture(posture,plotNose=true,IND=i,plot_handles=scatax);hold on
	% 	plotax=csp.plot3_posture(posture,plotNose=true,IND=i,plot_handles=plotax);
	% 	daspect([1,1,1])
	% 	drawnow;
	% end
end



%% save all data
xval_Dataset=csp2save.xval_Dataset_cleared;

save(fullfile(MasterDir,csp.currPlaneName4Directory+".mat"), ...
	"Both_score_ridge_cells","score_eigdim","xval_Dataset","ridgeparams");

disp("saving complete!")

toc(ticmaster);

%%

function Csp=saveExampleCellac(csp,scores_ridge_cell,numcells2save)
Csp=copy(csp);
maxscore_cell=max(scores_ridge_cell,[],1);
[~,bestcellind]=sort(maxscore_cell,"descend");
[~,bestridgeind]=max(scores_ridge_cell(:,bestcellind),[],1);
for g=1:numel(Csp.xval_Dataset)
	cellind=bestcellind(1:numcells2save);
	ridgeind=bestridgeind(1:numcells2save);
	IND={':',cellind,':',':'};
	Csp.xval_Dataset{g}.Tester.Responser=Csp.xval_Dataset{g}.Tester.Responser(IND{:});

	rec_ridge={};
	for i=1:numcells2save
		cellind=bestcellind(i);
		ridgeind=bestridgeind(i);
		IND={':',cellind,':',':'};
		rec_ridge{g}=Csp.xval_Dataset{g}.Tester.Reconstructed{ridgeind}(IND{:});
	end
	Csp.xval_Dataset{g}.Tester.Reconstructed=cat(2,rec_ridge{:});
end
Csp.organize_xval_Dataset(["Tester","Responser"],["Tester","Reconstructed"]);
Csp.xval_Dataset_cleared;
end







