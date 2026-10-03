%% assess direction selectivity and rotation, representative data
path2functions='/path/to/matlab_helpers';
path2allMiceDataDirlist="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2.mat";
path2BrainRegions="/path/to/data/BrainRegions.mat";
MasterDir="/path/to/results/example_neurons";	% path to store all data

addpath(genpath(path2functions));
%% direction selective cells
skewthresh=0.5;
pthresh=0.05;
numcells2plot=4;
linewidth=2;
colors=colororder();
for plane=18		% Animal_09
	tic
	currfolder=fullfile(MasterDir,CSP().PlaneName4Directory(plane));
	mkdir(currfolder);
	csp=CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist);
	csp.loadData(plane, ...
		PeriodCutoff4Derivatives_sec=0.1, ...
		triErrorThresh_mm=4);
	csp.setCurrentData("p8")
	%%
	csp.assess_direction_selectivity(divide_by_sum=false, ...
		pre_margin_in_sec=0);
	Csp=copy(csp);
	Csp.modify_trialTime(tstart_shift_sec=-1,tend_shift_sec=1);
	Skew=csp.Skew(:);
	[selectivity_sorted,selectivity_sorted_ind]=sort(csp.direction_selectivity.score,"descend");
	valid_sorted_ind=find((Skew(selectivity_sorted_ind)>=skewthresh)&...
		(csp.direction_selectivity.pval(selectivity_sorted_ind)<=0.05));
	%% plot selectivity score and pval
	figure;scatter(csp.direction_selectivity.score(csp.Skew>=skewthresh),csp.direction_selectivity.pval(csp.Skew>=skewthresh),5,"filled")
	%% 
	takeZscore=false;
	cat_8th_dir=true;
	plot_angle=true;
	cellinds=[4055,4364,4272];
	% cellinds=[4055,4364,4059];
	tckd=Csp.p8_tckd;
	if plot_angle
		ntl=nested_tiledlayout(5,1);
	else
		ntl=nested_tiledlayout(4,1);
	end
	colors=colororder();
	colors=colors([2,5,1],:);

	% plot angle
	if plot_angle
		ntl.nexttile();
		plot([0,360],[0,360],linewidth=2,Color=[0,0,0]);
		xticks([0:45:360])
		xticklabels(string([0:45:360])+"\circ");
		yticks([0:45:360])
		yticklabels(string([0:45:360])+"\circ");
		title("angle",FontSize=10)
	end

	ntl.nexttile();
	for i=1:numel(cellinds)
		cellind=cellinds(i);
		tck=csp.rot_tck(:,cellind,:);
		if takeZscore
			tck=zscore(tck,1,[1,3]);
		end
		tk=squeeze(tck);
		% plot(tk,color=0.7*[1,1,1]);hold on
		plotMeanSE(tk,1,Color=colors(i,:),ShadeAlpha=0.3);
		xticks(1:(size(tk,1)-1)/8:size(tk,1))
		xticklabels(string([0:45:360])+"\circ")
		hold on
		grid off
	end
	title("activity of rot",FontSize=15);
	% currxlim=gca().XLim;
	% d=(currxlim(2)-currxlim(1))/16;
	% xlim(currxlim+[-d,d])
	for i=1:numel(cellinds)
		ntl.nexttile();
		% cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=cellinds(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		if takeZscore
			tkd=zscore(tkd,1,[1,2,3]);
		end
		if cat_8th_dir
			tkd=cat(3,tkd(:,:,end),tkd);
			plotSideBySide(tkd,Linewidth=0.5,MeanColor=colors(i,:),ShadeAlpha=0.3,xticklabels=string([0:45:360])+"\circ",plotIndividualData=false);
		else
			plotSideBySide(tkd,Linewidth=0.5,MeanColor=colors(i,:),ShadeAlpha=0.3,xticklabels=string([45:45:360])+"\circ",plotIndividualData=false);
		end
		title(sprintf("cell # %d, p8",cellind),FontSize=15)
		grid off
	end



	%%
	tckd=Csp.p8_tckd;
	% tckd=tckd-mean(tckd(1:10,:,:,:),1);
	% ntl=nested_tiledlayout(3,4);
	ntl=nested_tiledlayout('flow');
	for i=1:10
		ntl.nexttile();
		% cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=selectivity_sorted_ind(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		plotSideBySide(tkd,Linewidth=0.5,ShadeAlpha=0.4,xticklabels=string([45:45:360]),xtickFontSize=15);
		title(sprintf("cell # %d, p8",cellind),FontSize=15)
		grid off

		ntl.nexttile();
		tck=csp.rot_tck(:,cellind,:);
		tk=squeeze(tck);
		plot(tk,color=0.7*[1,1,1]);hold on
		plotMeanSE(tk,1,Color=colors(1,:),ShadeAlpha=0.3);
		xticks(1:(size(tk,1)-1)/8:size(tk,1))
		xticklabels(string([0:45:360])+"\circ")
		grid off
		ylim([prctile(tk,1,"all"),prctile(tk,99,"all")])
		title(sprintf("cell # %d, rot",cellind),FontSize=15)
	end
	ntl.change_settings(Title=sprintf("%s, p8 & rot\nActivity of neurons with high direction selectivity",CSP().PlaneName(plane)),TitleSize=20)
	%% save
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, raw")).dsvg();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, raw")).png();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, raw")).fig();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, raw")).pdf();
	%% yaw
	tckd=Csp.p8_tckd;
	% tckd=tckd-mean(tckd(1:10,:,:,:),1);
	ntl=nested_tiledlayout(3,4);
	for i=1:6
		ntl.nexttile();
		% cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=selectivity_sorted_ind(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		plotSideBySide(tkd,Linewidth=0.5,ShadeAlpha=0.4,xticklabels=string([45:45:360]),xtickFontSize=15);
		title(sprintf("cell # %d, p8",cellind),FontSize=15)
		grid off

		ntl.nexttile();
		tck=csp.yaw_tck(:,cellind,:);
		tk=squeeze(tck);
		plot(tk,color=0.7*[1,1,1]);hold on
		plotMeanSE(tk,1,Color=colors(1,:),ShadeAlpha=0.3);
		tp=mean(csp.yaw_tpk(:,3,:),3);
		xticks(1:4:size(tp,1));
		xticklabels(compose("%.0f",rad2deg(tp(1:4:end)))+"\circ");
		grid off
		ylim([prctile(tk,1,"all"),prctile(tk,99,"all")])
		title(sprintf("cell # %d, yaw",cellind),FontSize=15)
	end
	ntl.change_settings(Title=sprintf("%s, p8 & yaw\nActivity of neurons with high direction selectivity",CSP().PlaneName(plane)),TitleSize=20)
	%% save
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, raw")).dsvg();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, raw")).png();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, raw")).fig();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, raw")).pdf();
	%% zvert
	tckd=Csp.p8_tckd;
	% tckd=tckd-mean(tckd(1:10,:,:,:),1);
	ntl=nested_tiledlayout(3,4);
	for i=1:6
		ntl.nexttile();
		% cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=selectivity_sorted_ind(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		plotSideBySide(tkd,Linewidth=0.5,ShadeAlpha=0.4,xticklabels=string([45:45:360]),xtickFontSize=15);
		title(sprintf("cell # %d, p8",cellind),FontSize=15)
		grid off

		ntl.nexttile();
		tck=csp.zvert_tck(:,cellind,:);
		tk=squeeze(tck);
		plot(tk,color=0.7*[1,1,1]);hold on
		plotMeanSE(tk,1,Color=colors(1,:),ShadeAlpha=0.3);
		tp=mean(csp.zvert_tpk(:,6,:),3);
		xticks(1:4:size(tp,1));
		xticklabels(compose("%.1f",tp(1:4:size(tp,1))));
		grid off
		ylim([prctile(tk,1,"all"),prctile(tk,99,"all")])
		title(sprintf("cell # %d, zvert",cellind),FontSize=15)
	end
ntl.change_settings(Title=sprintf("%s, p8 & zvert\nActivity of neurons with high direction selectivity",CSP().PlaneName(plane)),TitleSize=20)
%% save
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, raw")).dsvg();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, raw")).png();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, raw")).fig();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, raw")).pdf();
	
%% change parameters for direction selectivity assessment %%%%%%%%%%%%%
%%
	csp.assess_direction_selectivity(divide_by_sum=true, ...
		pre_margin_in_sec=1);
	Csp=copy(csp);
	Csp.modify_trialTime(tstart_shift_sec=-1,tend_shift_sec=1);
	Skew=csp.Skew(:);
	[selectivity_sorted,selectivity_sorted_ind]=sort(csp.direction_selectivity.score,"descend");
	valid_sorted_ind=find((Skew(selectivity_sorted_ind)>=skewthresh)&...
		(csp.direction_selectivity.pval(selectivity_sorted_ind)<=0.05));
	%% plot selectivity score and pval
	figure;scatter(csp.direction_selectivity.score(csp.Skew>=skewthresh),csp.direction_selectivity.pval(csp.Skew>=skewthresh),5,"filled")
	%% plot neural activity of selective cells
	tckd=Csp.p8_tckd;
	% tckd=tckd-mean(tckd(1:10,:,:,:),1);
	ntl=nested_tiledlayout(3,4);
	for i=1:6
		ntl.nexttile();
		cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=selectivity_sorted_ind(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		TDk=multiBlockCat(tkd,3,1);
		theta=[0:size(TDk,1)-1].*(2*pi/size(TDk,1))+2*pi/16;
		rho=TDk;
		theta(end+1)=theta(1);
		rho(end+1,:)=rho(1,:);
		polarplot(theta,rho,Color=0.7*[1,1,1]);hold on;
		x=[0:size(tkd,1):size(TDk,1)-1].*(2*pi/size(TDk,1))+2*pi/16;
		x=x+zeros(2,size(x,2));
		y=[zeros(1,size(x,2));ones(1,size(x,2))]*10;
		% polarplot(x,y,color="black")

		if i==1
			cp=plotCartInPolar();
		end
		rlim([prctile(TDk,10,"all"),prctile(TDk,99,"all")])
		cp.makeCartAxes(gca());
		[Mean,SE]=calculateMeanSE(rho,1,x=theta);
		[mean_x,mean_y]=cp.pol2cart(Mean.x,Mean.y);
		[SE_x,SE_y]=cp.pol2cart(SE.x,SE.y);
		plot(mean_x,mean_y,LineWidth=linewidth,Color=colors(1,:));
		patch(SE_x,SE_y,colors(1,:),facealpha=0.3,edgealpha=0);
		hold on
		scatter(0,0,20,"black","filled")
		cp.closeCartAxes();
		title(sprintf("cell # %d, p8",cellind),FontSize=15)

		% plot rot
		ntl.nexttile();
		tck=csp.rot_tck(:,cellind,:);
		tk=squeeze(tck);
		theta=[0:size(tk,1)-1].*(2*pi/size(tk,1));
		rho=tk;
		theta(end+1)=theta(1);
		rho(end+1,:)=rho(1,:);
		polarplot(theta,rho,color=0.7*[1,1,1]);
		if i==1
			rotpolar=plotCartInPolar();
		end
		rlim([prctile(rho,0,"all"),prctile(rho,99,"all")])
		rotpolar.makeCartAxes(gca());
		[Mean,SE]=calculateMeanSE(rho,1,x=theta);
		[mean_x,mean_y]=rotpolar.pol2cart(Mean.x,Mean.y);
		[SE_x,SE_y]=rotpolar.pol2cart(SE.x,SE.y);
		plot(mean_x,mean_y,LineWidth=linewidth,Color=colors(1,:));
		patch(SE_x,SE_y,colors(1,:),facealpha=0.3,edgealpha=0);
		hold on
		scatter(0,0,20,"black","filled")
		rotpolar.closeCartAxes();

		title(sprintf("cell # %d, rot",cellind),FontSize=15)
	end
	ntl.change_settings(Title=sprintf("%s, p8 & rot, skew>=%g, ANOVA p<=%g, pre-margin=1 sec, divide by sum = true\nActivity of neurons with high direction selectivity",CSP().PlaneName(plane),skewthresh,pthresh),TitleSize=20)
	%% save
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot polar, filtered")).dsvg();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot polar, filtered")).png();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot polar, filtered")).fig();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot polar, filtered")).pdf();
	%% rot
	tckd=Csp.p8_tckd;
	% tckd=tckd-mean(tckd(1:10,:,:,:),1);
	ntl=nested_tiledlayout(3,4);
	for i=1:6
		ntl.nexttile();
		cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=selectivity_sorted_ind(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		plotSideBySide(tkd,Linewidth=0.5,ShadeAlpha=0.4,xticklabels=string([45:45:360]),xtickFontSize=15);
		title(sprintf("cell # %d, p8",cellind),FontSize=15)
		grid off

		ntl.nexttile();
		tck=csp.rot_tck(:,cellind,:);
		tk=squeeze(tck);
		plot(tk,color=0.7*[1,1,1]);hold on
		plotMeanSE(tk,1,Color=colors(1,:),ShadeAlpha=0.3);
		xticks(1:(size(tk,1)-1)/8:size(tk,1))
		xticklabels(string([0:45:360]))
		grid off
		ylim([prctile(tk,1,"all"),prctile(tk,99,"all")])
		title(sprintf("cell # %d, rot",cellind),FontSize=15)
	end
	ntl.change_settings(Title=sprintf("%s, p8 & rot, skew>=%g, ANOVA p<=%g, pre-margin=1 sec, divide by sum = true\nActivity of neurons with high direction selectivity",CSP().PlaneName(plane),skewthresh,pthresh),TitleSize=20)
	%% save
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, filtered")).dsvg();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, filtered")).png();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, filtered")).fig();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and rot, filtered")).pdf();

	%% yaw
	tckd=Csp.p8_tckd;
	% tckd=tckd-mean(tckd(1:10,:,:,:),1);
	ntl=nested_tiledlayout(3,4);
	for i=1:6
		ntl.nexttile();
		% cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=selectivity_sorted_ind(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		plotSideBySide(tkd,Linewidth=0.5,ShadeAlpha=0.4,xticklabels=string([45:45:360]),xtickFontSize=15);
		title(sprintf("cell # %d, p8",cellind),FontSize=15)
		grid off

		ntl.nexttile();
		tck=csp.yaw_tck(:,cellind,:);
		tk=squeeze(tck);
		plot(tk,color=0.7*[1,1,1]);hold on
		plotMeanSE(tk,1,Color=colors(1,:),ShadeAlpha=0.3);
		tp=mean(csp.yaw_tpk(:,3,:),3);
		xticks(1:4:size(tp,1));
		xticklabels(compose("%.0f",rad2deg(tp(1:4:end)))+"\circ");
		grid off
		ylim([prctile(tk,1,"all"),prctile(tk,99,"all")])
		title(sprintf("cell # %d, yaw",cellind),FontSize=15)
	end
	ntl.change_settings(Title=sprintf("%s, p8 & yaw, skew>=%g, ANOVA p<=%g, pre-margin=1 sec, divide by sum = true\nActivity of neurons with high direction selectivity",CSP().PlaneName(plane),skewthresh,pthresh),TitleSize=20)
	%% save
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, filtered")).dsvg();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, filtered")).png();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, filtered")).fig();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and yaw, filtered")).pdf();
	%% zvert
	tckd=Csp.p8_tckd;
	% tckd=tckd-mean(tckd(1:10,:,:,:),1);
	ntl=nested_tiledlayout(3,4);
	for i=1:6
		ntl.nexttile();
		% cellind=selectivity_sorted_ind(valid_sorted_ind(i));
		cellind=selectivity_sorted_ind(i);
		tkd=squeeze(tckd(:,cellind,:,:));
		plotSideBySide(tkd,Linewidth=0.5,ShadeAlpha=0.4,xticklabels=string([45:45:360]),xtickFontSize=15);
		title(sprintf("cell # %d, p8",cellind),FontSize=15)
		grid off

		ntl.nexttile();
		tck=csp.zvert_tck(:,cellind,:);
		tk=squeeze(tck);
		plot(tk,color=0.7*[1,1,1]);hold on
		plotMeanSE(tk,1,Color=colors(1,:),ShadeAlpha=0.3);
		tp=mean(csp.zvert_tpk(:,6,:),3);
		xticks(1:4:size(tp,1));
		xticklabels(compose("%.1f",tp(1:4:size(tp,1))));
		grid off
		ylim([prctile(tk,1,"all"),prctile(tk,99,"all")])
		title(sprintf("cell # %d, zvert",cellind),FontSize=15)
	end
	ntl.change_settings(Title=sprintf("%s, p8 & zvert, skew>=%g, ANOVA p<=%g, pre-margin=1 sec, divide by sum = true\nActivity of neurons with high direction selectivity",CSP().PlaneName(plane),skewthresh,pthresh),TitleSize=20)
	%% save
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, filtered")).dsvg();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, filtered")).png();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, filtered")).fig();
	ImageSaver(gcf(),path=fullfile(currfolder,"p8 and zvert, filtered")).pdf();
	%%
	close all
	toc
	fprintf("%d / %d\n",plane, numel(CSP().PlaneName))
end







