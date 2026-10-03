%% Direction selectivity of single neurons during 8DT (Fig. 2f,g and Fig. 4b)
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 20);
path2functions='/path/to/matlab_helpers';
path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2.mat";
path2BrainRegions="/path/to/data/BrainRegions_with_RSC.mat";
MasterDir="/path/to/results/direction_selectivity";	% path to store all data

addpath(genpath(path2functions));
mkdir(MasterDir);


%% direction selectivity score, preferred direction and ANOVA p-value of every neuron
skewthresh=0.5;
pthresh=0.05;
divide_by_sum=true;	% mean responses normalized to sum to 1 across directions
pre_margin_in_sec=1;
dir_sel={};
skew_plane={};
RegLog_plane={};
ROI_plane={};
half_ROI_plane={};
rot_tc_plane={};
rot_tp_plane={};
for plane=1:numel(CSP().PlaneName)
	csp=CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
		dir_BrainRegions_used=path2BrainRegions);
	csp.loadData(plane, ...
		PeriodCutoff4Derivatives_sec=0.1, ...
		triErrorThresh_mm=4);
	csp.setCurrentData("p8")
	csp.assess_direction_selectivity(divide_by_sum=divide_by_sum, ...
		pre_margin_in_sec=pre_margin_in_sec);
	dir_sel{plane}=csp.direction_selectivity;
	skew_plane{plane}=csp.Skew;
	RegLog_plane{plane}=csp.Reg_inRegLog_reg;
	ROI_plane{plane}=csp.regROI;
	half_ROI_plane{plane}=csp.makeHalfCoord();

	% assess rot
	rot_tc_plane{plane}=mean(csp.rot_tck,3);

	% obtain rot selective direction
	rot_tp_plane{plane}=mean(csp.rot_tpk,3);

	fprintf("plane %d / %d\n",plane,numel(CSP().PlaneName))
end
%% save results to MasterDir
save(fullfile(MasterDir, "direction_selectivity_"+getTimestamp+".mat"), ...
    "dir_sel", ...
    "skew_plane", ...
    "RegLog_plane", ...
    "ROI_plane", ...
    "half_ROI_plane", ...
    "rot_tc_plane", ...
    "rot_tp_plane", ...
    "skewthresh", ...
    "pthresh", ...
    "divide_by_sum", ...
    "pre_margin_in_sec");
%% load data
% file written by the previous section
load("/path/to/results/direction_selectivity/direction_selectivity_<time>.mat")

%% Fig. 2f: preferred directions of the direction-selective neurons
dir_sel_thresh=0.1;
skew_thresh=0.5;
pval_thresh=0.05;
figure;
colors=colororder();
dircount_p_d=[];
for plane=1:numel(CSP().PlaneName)
	skewlog=skew_plane{plane}>skew_thresh;
	pvallog=dir_sel{plane}.pval(skewlog)<pval_thresh;
	score=dir_sel{plane}.score(skewlog);
	score_log=score>dir_sel_thresh;
	direction=dir_sel{plane}.direction(skewlog);
	validind=pvallog.*score_log;
	valdir=direction(logical(validind));
	temp=accumarray(valdir,1);
	dircount=zeros(8,1);dircount(1:numel(temp))=temp;
	dircount_p_d=cat(1,dircount_p_d,dircount(:)');
end
polarbar(sum(dircount_p_d,1),EdgeAlpha=1,FaceColor=colors(1,:),ThetaDir="counterclockwise",ThetaZeroLocation=90+45/2);
polarTicks8direction()
title(sprintf("# of neurons with direction selectivity>%g & p-val<%g, \nunder condition: skew>%g",dir_sel_thresh,pval_thresh,skew_thresh),FontSize=20)

%% Fig. 2g: percentage of direction-selective neurons in each area
dir_sel_thresh=0.1;
skew_thresh=0.5;
pval_thresh=0.05;
cell_count_thresh=50;
valcount_p_r=[];
allcount_p_r=[];
ntl=nested_tiledlayout();
colors=colororder();
for plane=1:numel(CSP().PlaneName)
	skewlog=skew_plane{plane}>skew_thresh;
	pvallog=dir_sel{plane}.pval(skewlog)<pval_thresh;
	score=dir_sel{plane}.score(skewlog);
	score_log=score>dir_sel_thresh;
	for reg=1:numel(RegLog_plane{plane})
		reglog=RegLog_plane{plane}{reg}(skewlog);
		curr_score_log=score_log(reglog);
		curr_pval_log=pvallog(reglog);
		curr_log=curr_score_log & curr_pval_log;
		valcount_p_r(plane,reg)=sum(curr_log);
		allcount_p_r(plane,reg)=numel(curr_log);
		if numel(curr_log)<cell_count_thresh
			valcount_p_r(plane,reg)=NaN;
			allcount_p_r(plane,reg)=NaN;
		end
	end
end
ntl.nexttile()

% -------- region order --------
desired_order = ["M2","M1","Whisker","UL","LL","Trunk","PPC","RSC"];
region_names_all = string(CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
	dir_BrainRegions_used=path2BrainRegions).RegionNames_type);

[tf, order_idx] = ismember(desired_order, region_names_all);
if ~all(tf)
    error("These requested regions were not found: %s", strjoin(desired_order(~tf), ", "));
end

% -------- ratio per plane x region --------
ratio_p_r = valcount_p_r ./ allcount_p_r;

% reorder columns
ratio_p_r    = ratio_p_r(:, order_idx);
region_names = region_names_all(order_idx);

% -------- mean ± SEM across planes --------
mu = mean(ratio_p_r, 1, "omitmissing");
sd = std(ratio_p_r, 0, 1, "omitmissing");
n  = sum(~isnan(ratio_p_r), 1);
se = sd ./ sqrt(n);

nReg = numel(region_names);

% Statistics (two-sided): direction_selectivity_stats.py

% -------- plot --------
bar(mu, EdgeAlpha=0, FaceColor=colors(1,:));
hold on

errorbar(1:nReg, mu, se, ...
    LineStyle="none", LineWidth=2, Color=colors(2,:))

text(1:nReg, mu, compose("%.2f", 100*mu) + "%", ...
    FontSize=13, ...
    HorizontalAlignment="center", ...
    VerticalAlignment="bottom")

xticks(1:nReg)
xticklabels(region_names)
xlim([0, nReg + 1])
ax = gca;
ax.YTickLabel = compose("%.0f%%", 100 * ax.YTick);

title_str = sprintf([ ...
    "Ratio of neurons with direction selectivity > %g & p-val < %g,\n" +...
    "under condition: skew > %g, neuron count threshold > %g"], ...
    dir_sel_thresh, pval_thresh, skew_thresh, cell_count_thresh);

ntl.change_settings( ...
    Xlabel="Brain Regions", ...
    Ylabel="Ratio of valid neurons", ...
    XlabelSize=20, ...
    YlabelSize=20, ...
    Title=title_str, ...
    TitleSize=20, ...
    LinkAxes='xy')

%% Fig. 4b: correlation between the 8DT and T-Rot responses of direction-selective neurons
skewthresh=0.5;
pthresh=0.05;
score_thresh=0.1;
corr_all=[];
for plane=1:numel(CSP().PlaneName)
	CD=dir_sel{plane}.CD;
	rot_tc=rot_tc_plane{plane};
	% interpolate p8 direction score	
	DC=CD';
	DCcat=[DC(end,:);DC];
	dlist=[0:45:360];
	rot_dlist=linspace(0,360,size(rot_tc,1));
	p8_time_cell=interp1(dlist,DCcat,rot_dlist);
	rot_time_cell=rot_tc;
	corr_cell=dot(zscore(p8_time_cell,0,1),zscore(rot_time_cell,0,1),1)./(size(p8_time_cell,1)-1);
	valid_log=(skew_plane{plane}(:)>skewthresh)&(dir_sel{plane}.pval(:)<pthresh)&(dir_sel{plane}.score(:)>score_thresh);
	corr_cell=corr_cell(valid_log);
	corr_all=[corr_all;corr_cell(:)];
end
figure;
histogram(corr_all,50,Normalization="probability");
xlabel("correlation between p8 and rot score")
ylabel("probability")
title(sprintf("Comparison between directin selectivity between p8 and rot" + ...
	"\nall neurons with score>%g, pval<%g, skew>%g",score_thresh,pthresh,skewthresh))
