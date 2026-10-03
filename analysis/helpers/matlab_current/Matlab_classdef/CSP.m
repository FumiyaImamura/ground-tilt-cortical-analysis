classdef CSP < matlab.mixin.Copyable
	% Stands for Cells, Stewart platform, Posture
    % This class uses info.NoAltplane.
    % Initialize the object by providing the info structure and a plane.

	% Example: 
	% for plane=1:length(info.NoAltplane.nameMouse_p)
	% Data=CSP(info=info);
	% Data.loadData(plane);
	% Data.setCurrentData("p8");
	% end
	properties (Constant)
		dir_BrainRegions ='/path/to/data/BrainRegions.mat'	% directory for brain regions
		dir_AllMiceDirectories="/path/to/data/AllMiceDirectories"	% directory for data location of each mice	%% Unused
		dir_AllMiceDataFile="/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat";
	end
    
    properties
        % Necessary input for initialization
        info    % Should be compatible with info from "/path/to/analysis_data/info_2025-02-05.mat"
		dirlist	% contains all the directories for all mice
        plane=[]   % Which plane (imaging session) to use
		currPlaneName
		currPlaneName4Directory

		% If these are specified, they override dir_BrainRegions and dir_AllMiceDataFile
		dir_BrainRegions_used	% directory for brain regions
		dir_AllMiceDataFile_used

		cellac
		org
		DLC_raw
		invDLC_raw

		cellac_mtype
		sp_mtype
		DLC_mtype
		invDLC_mtype

		SPparam_ind
		ind_SPparam
        
        % Raw data
        p8_tckd
        rot_tck
        tilt_tck
        yaw_tck
        zvert_tck

        p8_tpkd
        rot_tpk
        tilt_tpk
        yaw_tpk
        zvert_tpk

        p8_DLC
        rot_DLC
        tilt_DLC
        yaw_DLC
        zvert_DLC
        
        rec3d
		Skew
		regROI	% example of how they should be plotted: figure;imshow(CSP.AtlasImage);hold on;scatter(regROI(:,2),regROI(:,1));axis equal
				% Left top of the cortex image is the origin, and regROI(:,1) is the y axis (vertically down) and regROI(:,2) is the x axis (horizontally right) 

		AtlasImage
		RegionNames
		RegionNames_type
		contour_LeftRegion
		contour_LeftRegionLog
		contour_RightRegion
		contour_RightRegionLog
		contour_region

		% regions and roi for the current plane (hemisphere)
		Reg_LeftOrRight_reg
		Reg_contour_reg
		Reg_inRegLog_reg
		Reg_RegionNames

		%
		PlaneName
		PlaneName4Directory

		% plane ind and mouse name for each mouse
		plane_MouseIND
		mouseName_MouseIND

		SamplingRate_Hz
		PeriodCutoff4Derivatives_sec
		triErrorThresh
		DivideDerivWithTime

		% list of directories where the data are saved
		savedDirectories =struct
    end
    
    methods
        % Constructor: Note that the constructor name must match the class name.
		function obj = CSP(options)
            arguments
                options.?CSP  % Allow name–value pair arguments that match properties
            end
            % Loop through each provided property and assign it.
            props = fieldnames(options);
            for i = 1:numel(props)
                prop = props{i};
                obj.(prop) = options.(prop);
			end

			%%%%%%%%%% make info %%%%%%%%%%
			%%%%
			% datafolder=obj.dir_AllMiceDirectories;
			% datadir=dir(datafolder);
			% folderlist=string;
			% for i=3:length(datadir)
			% 	folderlist(i-2)=fullfile(datafolder,convertCharsToStrings(datadir(i).name));
			% end
			% allmicedir=[];
			% allmiceNames=string;
			% for i=1:length(folderlist)
			% 	temp=load(folderlist(i));
			% 	allmicedir{i}=temp.allDirectories;
			% 	allmiceNames(i)=temp.MouseName;
			% end
			%%%%

			%%%%
			if isempty(obj.dir_AllMiceDataFile_used)
				data=load(fullfile(obj.dir_AllMiceDataFile));
				obj.dir_AllMiceDataFile_used=obj.dir_AllMiceDataFile;
			else
				data=load(fullfile(obj.dir_AllMiceDataFile_used));
			end
			allmicedir=data.newallmicedir;
			allmiceNames=data.allmiceNames;
			%%%%

			temp=CellWithStruct2StructWithCell(allmicedir);
			obj.dirlist=structfun(@(a)string(a(:)),temp,UniformOutput=false);
			obj.PlaneName4Directory=allmiceNames(:);
			obj.PlaneName=regexprep(allmiceNames(:),"_"," ");

			info=[];
			info.NoAltplane=[];
			info.NoAltplane.nameFull_p=regexprep(allmiceNames(:),"_"," ");
			info.NoAltplane.nameMouse_p...
				=removeAllCharsAfterUnderScore(allmiceNames(:));
			info.NoAltplane.logRight_p=containsSubstring(allmiceNames(:),"right",CapSensitive=false);
			info.NoAltplane.logLeft_p=containsSubstring(allmiceNames(:),"left",CapSensitive=false);
			if any((info.NoAltplane.logRight_p+info.NoAltplane.logLeft_p)==0)
				error("Some planes cannot be assigned right or left")
			end
			obj.info=info;

			if isempty(obj.dir_BrainRegions_used)
				tempBrain = load(fullfile(obj.dir_BrainRegions));
				obj.dir_BrainRegions_used=obj.dir_BrainRegions;
			else
				tempBrain = load(fullfile(obj.dir_BrainRegions_used));
			end
			obj.AtlasImage = tempBrain.AtlasImage;
			obj.RegionNames = tempBrain.RegionNames;
			obj.RegionNames_type = tempBrain.RegionNames_type;
			obj.contour_LeftRegion = tempBrain.contour_LeftRegion;
			obj.contour_LeftRegionLog = tempBrain.contour_LeftRegionLog;
			obj.contour_RightRegion = tempBrain.contour_RightRegion;
			obj.contour_RightRegionLog = tempBrain.contour_RightRegionLog;
			obj.contour_region = tempBrain.contour_region;

			obj.SPparam_ind=["theta","phi","tau","X","Y","Z"];
			obj.ind_SPparam=@(b)arrayfun(@(a)find(obj.SPparam_ind==a),b);

			obj.generatePlaneIND4eachMouse();

			obj.calculatemidline();
			obj.makeHalfAtlas();

			function out=removeAllCharsAfterUnderScore(strArray)
				out = regexprep(strArray, '_.*', '');
			end
		end
	end
	%% load data
	properties
		RawQ_bp
	end
	methods
		function loadData(obj,plane,options)
			arguments
				obj
				plane=[]
				options.SamplingRate_Hz=10
				options.PeriodCutoff4Derivatives_sec=0.1
				options.smooth_posture_period_cutoff_sec=[]
				options.UseCachedData (1,1) logical =false	% true if loadData has been executed once, and you modified obj.cellac, obj.DLC_raw, or obj.org
				options.triErrorThresh_mm=[];
				options.Cached_cellac=[]
				options.Cached_org=[]
				options.Cached_DLC_raw=[]
				options.Cached_invDLC_raw=[]
				options.saveRawQ_bp=false
				options.DivideDerivWithTime=false
			end
			obj.SamplingRate_Hz=options.SamplingRate_Hz;
			obj.PeriodCutoff4Derivatives_sec=options.PeriodCutoff4Derivatives_sec;
			obj.triErrorThresh=options.triErrorThresh_mm;
			obj.DivideDerivWithTime=options.DivideDerivWithTime;
			if ~isempty(plane)
				if options.UseCachedData
					if obj.plane~=plane
						warning("plane and obj.plane is contradicting. Using obj.plane...")
						plane=obj.plane;
					end
				end
				obj.plane=plane;
			else
				if isempty(obj.plane)
					error("plane or obj.plane has to be specified")
				else
					plane=obj.plane;
				end
			end
			obj.currPlaneName=obj.PlaneName(obj.plane);
			obj.currPlaneName4Directory=obj.PlaneName4Directory(obj.plane);

			period_ms_cutoff=round(options.PeriodCutoff4Derivatives_sec*1000);

			if options.UseCachedData
				if isempty(obj.cellac)|isempty(obj.DLC_raw)|isempty(obj.org)
					error("Cellac,DLC_raw, or org is missing")
				end
			end

			% load all data
			if ~options.UseCachedData
				temp = load(fullfile(obj.dirlist.DLCdata_3DpointsOnly(plane)));	Q_bp = temp.Q_bp;
				temp = load(fullfile(obj.dirlist.DLCdata_info(plane)));			rec3d = temp.rec3d;
				temp = load(fullfile(obj.dirlist.SyncAndSP(plane)));              org   = temp.org;
				temp = load(fullfile(obj.dirlist.cellac(plane)));					cellac= temp.cellac;
				temp = load(fullfile(obj.dirlist.reg(plane)));					reg   = temp.reg;
				temp = load(fullfile(obj.dirlist.skewness(plane)));				skew  = temp.skew;
				
				% create function handle
				org.ind_mtypes=@(b)arrayfun(@(a)(find(org.mtypes==a)),b);
				rec3d.ind_bp=@(b)arrayfun(@(a)find(rec3d.bp_ind==a),b);
				rec3d.name_ind=rec3d.bp_ind;
				rec3d.ind_name=rec3d.ind_bp;
				rec3d.Dotname_ind=rec3d.jointList;
				rec3d.jointNames=join(rec3d.jointList,"_",2);
				org=fix_p8_org(org);
			else
				rec3d=obj.rec3d;
				skew=obj.Skew;
				reg=[];
				reg.roi=obj.regROI;
				if isempty(options.Cached_org)
					org=obj.org;
				else
					org=options.Cached_org;
				end
			end


			% organize cell activity
			if ~options.UseCachedData
				cellac		=NormalizeCellActivity(cellac);
				obj.cellac	=cellac;
			else
				if isempty(options.Cached_cellac)
					cellac=obj.cellac;
				else
					cellac=options.Cached_cellac;
				end
			end
			cellac_trial	=Sort_cellac_to_trials_AndDownsample(cellac,org,options.SamplingRate_Hz);
			cellac_mtype	=Sort_trials_to_mtypes(cellac_trial,org,3);

			% organize sp
			if ~options.UseCachedData
				org.sp_t	=PreprocessSPdata(org.sp_t);
				org			=fix_Rot_motion_org(org);
				org			=fix_p8_motion_org(org);

			end
			sp_trial		=Sort_sp_to_trials_AndDownsample(org,options.SamplingRate_Hz);
			sp_mtype		=Sort_trials_to_mtypes(sp_trial,org,3);
			if ~options.UseCachedData
				sp_mtype		=fixRotMotionSP(sp_mtype,org);
				sp_mtype		=fix_p8_MotionSP(sp_mtype,org);
			end


			% Organize DLC
			if options.saveRawQ_bp
				obj.RawQ_bp=Q_bp;
			end
			if ~options.UseCachedData
				DLC=[];
				DLC.tXYZbp					=Q_bp;
				DLC.triError				=cell2mat(rec3d.offset(:)');
				if ~isempty(options.triErrorThresh_mm)
					DLC=DLC_interpolateHighErrorDataPoints(DLC,options.triErrorThresh_mm);
				end
				[DLC.tANGbp,DLC.tDOTbp]		=DLC_make_Joint_Angle(Q_bp,rec3d);
				DLC							=DLC_calculateDerivatives(DLC,period_ms_cutoff,options.DivideDerivWithTime,org);
				DLC							=smooth_DLC(DLC,options.smooth_posture_period_cutoff_sec);
				obj.DLC_raw=DLC;

				% make invDLC %%%%%%%%%%%
				invDLC=[];
				invDLC.triError=DLC.triError;
				SP_T=SPparams2TransMat(org.sp_t);
				invDLC.tXYZbp=applyInverseSPTrans2TXYZbpUnsynced(SP_T,DLC.tXYZbp,org);
				[invDLC.tANGbp,invDLC.tDOTbp]=DLC_make_Joint_Angle(invDLC.tXYZbp,rec3d);
				invDLC							=DLC_calculateDerivatives(invDLC,period_ms_cutoff,options.DivideDerivWithTime,org);
				invDLC							=smooth_DLC(invDLC,options.smooth_posture_period_cutoff_sec);
				obj.invDLC_raw=invDLC;
				
			else
				if isempty(options.Cached_DLC_raw)
					DLC=obj.DLC_raw;
					invDLC=obj.invDLC_raw;
				else
					DLC=options.Cached_DLC_raw;
					invDLC=options.Cached_invDLC_raw;
				end
			end
			DLC_trial					=Sort_DLC_to_trials_AndDownsample(DLC,org,options.SamplingRate_Hz);
			DLC_mtype					=Sort_DLC_trials_to_mtypes(DLC_trial,org);
			
			invDLC_trial					=Sort_DLC_to_trials_AndDownsample(invDLC,org,options.SamplingRate_Hz);
			invDLC_mtype					=Sort_DLC_trials_to_mtypes(invDLC_trial,org);
			
			obj.cellac_mtype	=cellac_mtype;
			obj.sp_mtype		=sp_mtype;
			obj.DLC_mtype		=DLC_mtype;
			obj.invDLC_mtype	=invDLC_mtype;

			if assessNaN(obj.cellac_mtype);error("cellac has nan");end
			if assessNaN(obj.sp_mtype);error("sp has nan");end
			if assessNaN(obj.DLC_mtype);error("DLC had nan");end
			
			[obj.p8_tckd,obj.p8_tpkd,obj.p8_DLC]		=make_p8_data();
			[obj.rot_tck,obj.rot_tpk,obj.rot_DLC]		=make_mtype_data("rot");
			[obj.tilt_tck,obj.tilt_tpk,obj.tilt_DLC]	=make_mtype_data("tilt");
			[obj.yaw_tck,obj.yaw_tpk,obj.yaw_DLC]		=make_mtype_data("yaw");
			[obj.zvert_tck,obj.zvert_tpk,obj.zvert_DLC]	=make_mtype_data("zvert");

			if ~options.UseCachedData
				obj.rec3d=rec3d;
				obj.Skew=skew;
				obj.regROI=reg.roi;
				obj.org=org;
			end

			% assess hemisphere
			obj.assessInsideBrainRegions();

			obj.calculatemidline();
			obj.makeHalfAtlas();
			obj.makeHalfCoord();

			function normalized_cellac=NormalizeCellActivity(cellac)
				normalized_cellac=(cellac-mean(cellac,1))./std(cellac,0,1);
			end
			function org=fix_p8_org(org)
				% check if motionlog is the same size as the t_strial
				if numel(org.motionLog)~=numel(org.t_strial)
					error("motionLog is supposed to have the same number of elements with t_strial")
				end
				interTrialSec=1;
				interTrialMilSec=interTrialSec*1000;
				ValidMotionInds=find(containsSubstring(org.mtypes,"p8"));
				for trial=1:length(org.motionLog)-1
					if ismember(org.motionLog(trial),ValidMotionInds)	% assess if current trial is p8 or not
						if ismember(org.motionLog(trial+1),ValidMotionInds)	% assess if next trial is p8 or not
							currTrialTheta=org.sp_t(org.t_strial(trial):org.t_etrial(trial),1);
							currTrialThetaEx=org.sp_t(org.t_strial(trial):org.t_strial(trial+1)-1,1);
							currTrialThetaEx=mode(currTrialTheta,1)+zeros(size(currTrialThetaEx));
							org.sp_t(org.t_strial(trial):org.t_strial(trial+1)-1,1)=currTrialThetaEx;
						end
					end
				end
			end
			function sp_mtype=fixRotMotionSP(sp_mtype,org)
				ind=org.ind_mtypes("rot");
				tpk=sp_mtype{ind}(:,1,:);

				for i=1:size(tpk,3)
					tpk(:,:,i)=force_monotonic_increase(tpk(:,:,i));
				end

				% % Assume tpk is [T, P, N] (time x parameters x trials)
				% T = size(tpk,1);
				% N = size(tpk,3);
				% 
				% % 1. Determine where all parameters are increasing between consecutive time points.
				% inc = squeeze(all(diff(tpk, 1, 1) > 0, 2));  % Now inc is [T-1, N]
				% 
				% % 2. Append false at the beginning and end for each trial, then compute diff to mark patch boundaries.
				% inc_ext = [false(1,N); inc; false(1,N)];  % now size: [T+1, N]
				% d = diff(double(inc_ext));                % d==1 marks start, d==-1 marks end; size: [T, N]
				% 
				% % 3. For each trial, find start and end indices of all increasing patches.
				% startCells = arrayfun(@(j) find(d(:,j)==1), 1:N, 'UniformOutput', false);
				% endCells   = arrayfun(@(j) find(d(:,j)==-1)-1, 1:N, 'UniformOutput', false);
				% 
				% % For each trial, choose the longest patch.
				% [~, idxLongest] = cellfun(@(s,e) max(e-s+1), startCells, endCells);
				% startTime = cellfun(@(s,i) s(i), startCells, num2cell(idxLongest));
				% endTime   = cellfun(@(e,i) e(i)+1, endCells, num2cell(idxLongest)); % add one because diff is one less than tpk
				% 
				% % 4. Create a logical mask that is true only for the longest patch in each trial.
				% tvec = (1:T)';  % column vector of time indices
				% mask = (tvec >= startTime) & (tvec <= endTime);  % uses implicit expansion; size: [T, N]
				% mask = reshape(mask, T, 1, N);                   % reshape to 3D: [T, 1, N]
				% mask = repmat(mask, [1, size(tpk,2), 1]);         % replicate across the parameter dimension
				% 
				% % 5. Set values outside the longest monotonic increasing patch to NaN.
				% tpk(~mask) = NaN;
				% tpk=fillmissing(tpk,"previous",1);
				% tpk=fillmissing(tpk,"next",1);
				% 
				% if any(isnan(tpk),"all"); error("There are NaN values in tpk");end
				% if any(diff(tpk,1,1)<0)
				% 	warning("There still remains timepoints where phi is not monotonically increasing");
				% end
				sp_mtype{ind}(:,1,:)=tpk;
			end
			function sp_mtype=fix_p8_MotionSP(sp_mtype,org)
				% make it so that theta is constant across the trial
				inds=containsSubstring(org.mtypes,"p8");
				inds=find(inds);
				for ind=inds(:)'
					tpk=sp_mtype{ind}(:,1,:,:);
					tpk=mode(tpk,1)+zeros(size(tpk));
					sp_mtype{ind}(:,1,:,:)=tpk;
				end
			end
			function org=fix_Rot_motion_org(org)
				for i=1:numel(org.t_strial)
					if org.motionLog(i)==1
						if i~=numel(org.t_strial)
							theta=org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1);
						else
							theta=org.sp_t(org.t_strial(i):org.t_etrial(i),1);
						end
						theta=wrapTo360(rad2deg(theta));
						theta=force_monotonic_increase(theta, min=0, max=360);
						theta=wrapTo360(theta-45/2)+45/2;
						theta=deg2rad(theta);
						if i~=numel(org.t_strial)
							org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1)=theta;
						else
							org.sp_t(org.t_strial(i):org.t_etrial(i),1)=theta;
						end
					end
				end
				for i=1:numel(org.t_strial)-1
					if isnan(org.motionLog(i)) && (org.motionLog(i+1)==1)
						if i~=numel(org.t_strial)
							theta=org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1);
						else
							theta=org.sp_t(org.t_strial(i):org.t_etrial(i),1);
						end
						theta=wrapTo360(rad2deg(theta));
						theta=force_monotonic_increase(theta, min=0, max=360);
						theta=wrapTo360(theta-45/2)+45/2;
						theta=deg2rad(theta);
						if i~=numel(org.t_strial)
							org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1)=theta;
						else
							org.sp_t(org.t_strial(i):org.t_etrial(i),1)=theta;
						end
					end
				end
			end
			function org=fix_p8_motion_org(org)
				for i=1:numel(org.t_strial)
					if org.motionLog(i)>=5
						if i~=numel(org.t_strial)
							theta=org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1);
						else
							theta=org.sp_t(org.t_strial(i):org.t_etrial(i),1);
						end
						theta=median(theta(:))+zeros(size(theta));
						if i~=numel(org.t_strial)
							org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1)=theta;
						else
							org.sp_t(org.t_strial(i):org.t_etrial(i),1)=theta;
						end
					end
				end
			end
			function cellac_trial=Sort_cellac_to_trials_AndDownsample(cellac,org,SamplingRate_Hz)
				cellac_trial=[];
				downsampled_t=make_downsampled_t_4eachTrials(org,SamplingRate_Hz);
				% for i=1:length(org.t_strial)
				% 	f2p_ind=org.f2p_t(downsampled_t{i});
				% 	cellac_trial{i}=interp1([1:size(cellac,1)],cellac,f2p_ind);
				% end
				vec_t=cell2mat(downsampled_t(:)')';
				f2p_ind=org.f2p_t(vec_t);
				vec_cellac=interp1([1:size(cellac,1)],cellac,f2p_ind);
				t_length=cellfun(@length,downsampled_t,UniformOutput=true);
				cellac_trial=mat2cell(vec_cellac,t_length,size(cellac,2));
			end
			function organized_data=Sort_trials_to_mtypes(trial_sorted_cellArray,org,dim2concat)
				data=trial_sorted_cellArray(:)';
				organized_data=cell(1,numel(org.mtypes));
				for i=1:numel(organized_data)
					ind=find(org.motionLog==i);
					celldata=data(ind);
					organized_data{i}=cat(dim2concat,celldata{:});
				end
				% for i=1:numel(data)
				% 	if ~isnan(org.motionLog(i))
					% 	organized_data{org.motionLog(i)}=cat(dim2concat,organized_data{org.motionLog(i)},data{i});
				% 	end
				% end
			end
			function sp_trial=Sort_sp_to_trials_AndDownsample(org,SamplingRate_Hz)
				sp_trial=[];
				downsampled_t=make_downsampled_t_4eachTrials(org,SamplingRate_Hz);
				for i=1:length(org.t_strial)
					sp_trial{i}=org.sp_t(downsampled_t{i},:);
				end
			end
			function DLC=DLC_interpolateHighErrorDataPoints(DLC,triErrorThresh)
				Qind1=cellfun(@(a)size(a,1),DLC.tXYZbp);
				Qind2=cellfun(@(a)size(a,2),DLC.tXYZbp);
				Q=cell2mat(DLC.tXYZbp(:));
				triError=DLC.triError(:);
				lowErrorlog=triError<=triErrorThresh;
				% make sure that the endpoints of each body part do not get
				% mixed up in interpolation
				for i=1:length(DLC.tXYZbp)
					ind=(i-1)*Qind1(1)+1:i*Qind1(1);
					a=Q(ind,:);
					lowE=lowErrorlog(ind);
					if lowE(1)==0
						a(1,:)=a(find(lowE,1,"first"),:);
						lowE(1)=1;
					end
					if lowE(end)==0
						a(end,:)=a(find(lowE,1,"last"),:);
						lowE(end)=1;
					end
					lowErrorlog(ind)=lowE;
					Q(ind,:)=a;
				end
				x=find(lowErrorlog);
				y=Q(lowErrorlog,:);
				xq=[1:size(Q,1)]';
				interpy=interp1(x,y,xq);
				interpy=fillmissing(interpy,"previous",1);
				interpy=fillmissing(interpy,"next",1);
				newQ=mat2cell(interpy,Qind1,size(interpy,2));
				newQ=reshape(newQ,size(DLC.tXYZbp));
				DLC.tXYZbp=newQ;
			end
			function [DLCang,DLCdot]=DLC_make_Joint_Angle(Q_bp,rec3d)
				DLCang=[];
				DLCdot=[];
				for i=1:size(rec3d.jointList,1)
					inds=rec3d.ind_bp(rec3d.jointList(i,:));
					Q=Q_bp(inds);
					a=Q{2}-Q{1};	a=a./sqrt(dot(a,a,2));
					b=Q{3}-Q{1};	b=b./sqrt(dot(b,b,2));
					ang=acos(dot(a,b,2));
					DLCang(:,i)=ang;
					DLCdot(:,i)=dot(a,b,2);
				end
			end
			function DLC=DLC_calculateDerivatives(DLC,period_ms_cutoff,DivideDerivWithTime,org)
				DLC.tDXYZbp=[];
				DLC.tD2XYZbp=[];
				DLCsamplingRate=1000*numel(org.t_fdlc)./(org.t_fdlc(end)-org.t_fdlc(1));
				for i=1:length(DLC.tXYZbp)
					[DLC.tDXYZbp{i},DLC.tD2XYZbp{i},filtered_data]=calculateDerivatives_along1stdim(DLC.tXYZbp{i},period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate);
					%%%%%%%% plot filtered data
					if false
						figure;plot(DLC.tXYZbp{i}(:,1),LineWidth=2);hold on;plot(filtered_data(:,1),LineWidth=2);
						xticks(0:1000:size(filtered_data,1))
						xlabel("time (ms)")
						legend("original data","filtered data")
						title("all the periods shorter than 100 ms is removed")
						ChangeFigureSettings(TitleSize=20,XLabelSize=15,LegendSize=15)
					end
				end
				[DLC.tDANGbp,DLC.tD2ANGbp,filtered_data]=calculateDerivatives_along1stdim(DLC.tANGbp,period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate);
				[DLC.tDDOTbp,DLC.tD2DOTbp,filtered_data]=calculateDerivatives_along1stdim(DLC.tDOTbp,period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate);
			end
			function smoothed_DLC=smooth_DLC(DLC,period_cutoff_sec)
				if ~isempty(period_cutoff_sec)
					smoothed_DLC=DLC;
					for i=1:length(DLC.tXYZbp)
						smoothed_DLC.tXYZbp{i}=lowpass_filter(DLC.tXYZbp{i},300,1/period_cutoff_sec);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sampling rate is 300, not 1000	2025/06/29
					end
					smoothed_DLC.tANGbp=lowpass_filter(DLC.tANGbp,300,1/period_cutoff_sec);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sampling rate is 300, not 1000	2025/06/29
					smoothed_DLC.tDOTbp=lowpass_filter(DLC.tDOTbp,300,1/period_cutoff_sec);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sampling rate is 300, not 1000	2025/06/29
				else
					smoothed_DLC=DLC;
				end
			end
			function [ddata,d2data,filtered_data]=calculateDerivatives_along1stdim(data,period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate)
				cutoff=floor(size(data,1)/period_ms_cutoff)+1;
				% filtered_data=lowpass_filter_mat_along1stdim(data,cutoff,Display_WarningIfNonReal=false);
				filtered_data=lowpass_filter(data,300,1000/period_ms_cutoff);	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sampling rate is 300, not 1000	2025/06/29
				ddata=diff(filtered_data,1,1);
				if DivideDerivWithTime
					ddata=ddata./(1/DLCsamplingRate);
				end

				% changed acceleration derivation 2025/06/29 %%%%%%%%%%%%%%
				% d2data=diff(filtered_data,2,1);
				sampling_rate=DLCsamplingRate;
				diff_sec=0.1;
				diff_time_points=round(sampling_rate*diff_sec);
				d2data=ddata(diff_time_points:end,:)-ddata(1:end-diff_time_points+1,:);
				d2data(size(d2data,1)+1:size(ddata,1),:)=NaN;
				d2data=fillmissing(d2data,"previous",1);
				if DivideDerivWithTime
					d2data=d2data./(diff_sec);
				end
			end
			function DLC_trial=Sort_DLC_to_trials_AndDownsample(DLC,org,SamplingRate_Hz)
				DLC_trial=[];
				downsampled_t=make_downsampled_t_4eachTrials(org,SamplingRate_Hz);
				DLC_trial=applyFunctionLayers(DLC,@(a)DownSample_and_Sort_mat_2_Trials(a,downsampled_t));
				function matcell=DownSample_and_Sort_mat_2_Trials(mat,downsampled_t)
					matcell=[];
					for i=1:length(downsampled_t)
						fdlc_ind=org.fdlc_t(downsampled_t{i});
						fdlc_ind=round(fdlc_ind);
						newmat=mat(fdlc_ind,:);
						matcell{i}=newmat;
					end
					% vec_t=cell2mat(downsampled_t(:)')';
					% fdlc_ind=org.fdlc_t(vec_t);
					% vec_mat=interp1([1:size(mat,1)],mat,fdlc_ind);
					% t_length=cellfun(@length,downsampled_t,UniformOutput=true);
					% matcell=mat2cell(vec_mat,t_length,size(mat,2));
				end
			end
			function DLC_mtype=Sort_DLC_trials_to_mtypes(DLC_trial,org)
				flds=string(fields(DLC_trial));
				xyzflds=flds(containsSubstring(flds,"xyz",CapSensitive=false));
				scalarflds=flds(containsSubstring(flds,"ang",CapSensitive=false) ...
					| containsSubstring(flds,"dot",CapSensitive=false) ...
					| containsSubstring(flds,"triError",CapSensitive=false));
				DLC_mtype=sort_to_mtypes(DLC_trial);
				function DLC_mtype=sort_to_mtypes(DLC_trial)
					sortdata=[];
					f=@(a)Sort_trials_to_mtypes(a,org,3);
					for i=xyzflds(:)'
						sortdata.(i)=cellfun(f,DLC_trial.(i),UniformOutput=false);
					end
					for i=scalarflds(:)'
						sortdata.(i)=f(DLC_trial.(i));
					end
					DLC_mtype=[];
					for i=1:length(org.mtypes)
						DLC_mtype{i}=[];
						for j=xyzflds(:)'
							DLC_mtype{i}.(j)=cellfun(@(a)a{i},sortdata.(j),UniformOutput=false);
						end
						for j=scalarflds(:)'
							DLC_mtype{i}.(j)=sortdata.(j){i};
						end
					end
				end
			end
			function [p8_tckd,p8_tpkd,p8_DLC]=make_p8_data()
				ind=find(containsSubstring(org.mtypes,"p8"));
				p8_tckd=make_cellac();
				p8_tpkd=make_sp();
				p8_DLC=make_DLC();

				function p8_tckd=make_cellac()
					data=obj.cellac_mtype(ind);
					sizearray=cellfun(@(a)size(a,3),data,UniformOutput=true);
					minsize=min(sizearray);
					aligneddata=cellfun(@(a)a(:,:,1:minsize),data,UniformOutput=false);
					p8_tckd=cat(4,aligneddata{:});
				end
				function p8_tpkd=make_sp()
					data=obj.sp_mtype(ind);
					sizearray=cellfun(@(a)size(a,3),data,UniformOutput=true);
					minsize=min(sizearray);
					aligneddata=cellfun(@(a)a(:,:,1:minsize),data,UniformOutput=false);
					p8_tpkd=cat(4,aligneddata{:});
				end
				function p8_DLC=make_DLC()
					data=obj.DLC_mtype(ind);
					sizearray=cellfun(@(a)size(a,3),obj.cellac_mtype(ind),UniformOutput=true);
					minsize=min(sizearray);
					flds=string(fields(data{1}));
					xyzflds=flds(containsSubstring(flds,"xyz",CapSensitive=false));
					scalarflds=flds(containsSubstring(flds,"ang",CapSensitive=false) ...
						| containsSubstring(flds,"dot",CapSensitive=false) ...
						| containsSubstring(flds,"triError",CapSensitive=false));
					data=make_data_XYZ_numericArray(data);
					p8_DLC=struct;
					for j=flds(:)'
						p8_DLC.(j)=[];
						for i=1:length(data)
							p8_DLC.(j)=cat(4,p8_DLC.(j),data{i}.(j)(:,:,1:minsize));
						end
					end
					p8_DLC=make_data_XYZ_cellArray(p8_DLC);

					function data=make_data_XYZ_numericArray(data)
						for i=1:length(data)
							for j=xyzflds(:)'
								temp=data{i}.(j);
								data{i}.(j)=cell2mat(temp(:)');
							end
						end
					end
					function data=make_data_XYZ_cellArray(data)
						for j=xyzflds(:)'
							a=[];
							for k=1:size(data.(j),2)/3
								IND=createIND(data.(j),2,3*(k-1)+1:3*k);
								a{k}=data.(j)(IND{:});
							end
							data.(j)=a;
						end
					end
				end
			end
			function [mtype_tck,mtype_tpk,mtype_DLC]=make_mtype_data(str)
				ind=org.ind_mtypes(str);
				mtype_tck=obj.cellac_mtype{ind};
				mtype_tpk=obj.sp_mtype{ind};
				mtype_DLC=obj.DLC_mtype{ind};
			end
			function downsampled_t=make_downsampled_t_4eachTrials(org,SamplingRate_Hz)
				downsampled_t=[];
				for i=1:length(org.t_strial)
					t=org.t_strial(i):org.t_etrial(i);
					t_downsampled=interp1([1:length(t)],t,[1:1000/SamplingRate_Hz:length(t)]);
					t_downsampled=round(t_downsampled);
					downsampled_t{i}=t_downsampled;
				end
			end
		end
	end
	%% load data for continuous time
	properties
		timestamp
		trials
		cellac_continuous
		sp_continuous
		DLC_countinuous
		invDLC_countinuous
	end
	methods
		function loadData_continuous_time(obj,plane,options)
			arguments
				obj
				plane=[]
				options.SamplingRate_Hz=10
				options.PeriodCutoff4Derivatives_sec=0.1
				options.smooth_posture_period_cutoff_sec=[]
				options.UseCachedData (1,1) logical =false	% true if loadData has been executed once, and you modified obj.cellac, obj.DLC_raw, or obj.org
				options.triErrorThresh_mm=[];
				options.Cached_cellac=[]
				options.Cached_org=[]
				options.Cached_DLC_raw=[]
				options.Cached_invDLC_raw=[]
				options.DivideDerivWithTime=false
			end
			if options.UseCachedData
				options.SamplingRate_Hz=obj.SamplingRate_Hz;
				options.PeriodCutoff4Derivatives_sec=obj.PeriodCutoff4Derivatives_sec;
				options.triErrorThresh_mm=obj.triErrorThresh;
				options.DivideDerivWithTime=obj.DivideDerivWithTime;
			else
				obj.SamplingRate_Hz=options.SamplingRate_Hz;
				obj.PeriodCutoff4Derivatives_sec=options.PeriodCutoff4Derivatives_sec;
				obj.triErrorThresh=options.triErrorThresh_mm;
				obj.DivideDerivWithTime=options.DivideDerivWithTime;
			end
			if ~isempty(plane)
				if options.UseCachedData
					if obj.plane~=plane
						warning("plane and obj.plane is contradicting. Using obj.plane...")
						plane=obj.plane;
					end
				end
				obj.plane=plane;
			else
				if isempty(obj.plane)
					error("plane or obj.plane has to be specified")
				else
					plane=obj.plane;
				end
			end
			obj.currPlaneName=obj.PlaneName(obj.plane);
			obj.currPlaneName4Directory=obj.PlaneName4Directory(obj.plane);

			period_ms_cutoff=round(options.PeriodCutoff4Derivatives_sec*1000);

			if options.UseCachedData
				if isempty(obj.cellac)|isempty(obj.DLC_raw)|isempty(obj.org)
					error("Cellac,DLC_raw, or org is missing")
				end
			end

			% load all data
			if ~options.UseCachedData
				temp = load(fullfile(obj.dirlist.DLCdata_3DpointsOnly(plane)));	Q_bp = temp.Q_bp;
				temp = load(fullfile(obj.dirlist.DLCdata_info(plane)));			rec3d = temp.rec3d;
				temp = load(fullfile(obj.dirlist.SyncAndSP(plane)));              org   = temp.org;
				temp = load(fullfile(obj.dirlist.cellac(plane)));					cellac= temp.cellac;
				temp = load(fullfile(obj.dirlist.reg(plane)));					reg   = temp.reg;
				temp = load(fullfile(obj.dirlist.skewness(plane)));				skew  = temp.skew;
				
				% create function handle
				org.ind_mtypes=@(b)arrayfun(@(a)(find(org.mtypes==a)),b);
				rec3d.ind_bp=@(b)arrayfun(@(a)find(rec3d.bp_ind==a),b);
				rec3d.name_ind=rec3d.bp_ind;
				rec3d.ind_name=rec3d.ind_bp;
				rec3d.Dotname_ind=rec3d.jointList;
				rec3d.jointNames=join(rec3d.jointList,"_",2);
				% org=fixROT(org);
			else
				rec3d=obj.rec3d;
				skew=obj.Skew;
				reg=[];
				reg.roi=obj.regROI;
				if isempty(options.Cached_org)
					org=obj.org;
				else
					org=options.Cached_org;
				end
			end

			% make timestamps for downsampling
			timestamp=[1:1000/options.SamplingRate_Hz:length(org.f2p_t)];
			timestamp=round(timestamp);	% --> integer timpoint
			obj.timestamp=timestamp;


			% organize cell activity
			if ~options.UseCachedData
				cellac		=NormalizeCellActivity(cellac);
				obj.cellac	=cellac;
			else
				if isempty(options.Cached_cellac)
					cellac=obj.cellac;
				else
					cellac=options.Cached_cellac;
				end
			end
			cellac_downsampled=Downsample_cellac(cellac,timestamp,org);
			% cellac_trial	=Sort_cellac_to_trials_AndDownsample(cellac,org,options.SamplingRate_Hz);
			% cellac_mtype	=Sort_trials_to_mtypes(cellac_trial,org,3);

			% organize sp
			if ~options.UseCachedData
				prevorg=org;
				org.sp_t	=PreprocessSPdata_Accurate(org);
				raw_sp_t=org.sp_t;
				org=fix_p8(org,prevorg);
				org			=fix_Rot_motion_org(org);
			end
			sp_downsampled=org.sp_t(timestamp,:);
			% sp_trial% 		=Sort_sp_to_trials_AndDownsample(org,options.SamplingRate_Hz);
			% sp_mtype% 		=Sort_trials_to_mtypes(sp_trial,org,3);
			% sp_mtype% 		=fixRotMotionSP(sp_mtype,org);


			% Organize DLC
			if ~options.UseCachedData
				DLC=[];
				DLC.tXYZbp					=Q_bp;
				DLC.triError				=cell2mat(rec3d.offset(:)');
				if ~isempty(options.triErrorThresh_mm)
					DLC=DLC_interpolateHighErrorDataPoints(DLC,options.triErrorThresh_mm);
				end
				[DLC.tANGbp,DLC.tDOTbp]		=DLC_make_Joint_Angle(Q_bp,rec3d);
				DLC							=DLC_calculateDerivatives(DLC,period_ms_cutoff,options.DivideDerivWithTime,org);
				DLC							=smooth_DLC(DLC,options.smooth_posture_period_cutoff_sec);
				obj.DLC_raw=DLC;

				% make invDLC %%%%%%%%%%%
				invDLC=[];
				invDLC.triError=DLC.triError;
				SP_T=SPparams2TransMat(raw_sp_t);
				invDLC.tXYZbp=applyInverseSPTrans2TXYZbpUnsynced(SP_T,DLC.tXYZbp,org);
				[invDLC.tANGbp,invDLC.tDOTbp]=DLC_make_Joint_Angle(invDLC.tXYZbp,rec3d);
				invDLC							=DLC_calculateDerivatives(invDLC,period_ms_cutoff,options.DivideDerivWithTime,org);
				invDLC							=smooth_DLC(invDLC,options.smooth_posture_period_cutoff_sec);
				obj.invDLC_raw=invDLC;
			else
				if isempty(options.Cached_DLC_raw)
					DLC=obj.DLC_raw;
					invDLC=obj.invDLC_raw;
				else
					DLC=options.Cached_DLC_raw;
					invDLC=options.Cached_invDLC_raw;
				end
			end
			DLC_downsampled=Downsample_DLC(DLC,timestamp,org);
			invDLC_downsampled=Downsample_DLC(invDLC,timestamp,org);
			% DLC_trial					=Sort_DLC_to_trials_AndDownsample(DLC,org,options.SamplingRate_Hz);
			% DLC_mtype					=Sort_DLC_trials_to_mtypes(DLC_trial,org);

			trials=Downsample_trial_border(org,timestamp);


			obj.trials				=trials;
			obj.cellac_continuous	=cellac_downsampled;
			obj.sp_continuous		=sp_downsampled;
			obj.DLC_countinuous		=DLC_downsampled;
			obj.invDLC_countinuous	=invDLC_downsampled;
			

			
			% obj.cellac_mtype	=cellac_mtype;
			% obj.sp_mtype		=sp_mtype;
			% obj.DLC_mtype		=DLC_mtype;
			% 
			% if assessNaN(obj.cellac_mtype);error("cellac has nan");end
			% if assessNaN(obj.sp_mtype);error("sp has nan");end
			% if assessNaN(obj.DLC_mtype);error("DLC had nan");end
			% 
			% [obj.p8_tckd,obj.p8_tpkd,obj.p8_DLC]		=make_p8_data();
			% [obj.rot_tck,obj.rot_tpk,obj.rot_DLC]		=make_mtype_data("rot");
			% [obj.tilt_tck,obj.tilt_tpk,obj.tilt_DLC]	=make_mtype_data("tilt");
			% [obj.yaw_tck,obj.yaw_tpk,obj.yaw_DLC]		=make_mtype_data("yaw");
			% [obj.zvert_tck,obj.zvert_tpk,obj.zvert_DLC]	=make_mtype_data("zvert");

			if ~options.UseCachedData
				obj.rec3d=rec3d;
				obj.Skew=skew;
				obj.regROI=reg.roi;
				obj.org=org;
			end

			% assess hemisphere
			obj.assessInsideBrainRegions();

			obj.calculatemidline();
			obj.makeHalfAtlas();
			obj.makeHalfCoord();

			function normalized_cellac=NormalizeCellActivity(cellac)
				normalized_cellac=(cellac-mean(cellac,1))./std(cellac,0,1);
			end
			function org=fix_p8(org,prevorg)
				theta=prevorg.sp_t(:,1);
				newtheta=org.sp_t(:,1);
				costheta=org.cossp_t(:,1);
				sintheta=org.sinsp_t(:,1);
				ind=find(containsSubstring(org.mtypes,"p8"));
				for i=1:length(org.t_strial)
					if ~ismember(org.motionLog(i),ind)
						continue;
					end
					if i<length(org.t_strial)
						IND=org.t_strial(i):org.t_strial(i+1);
					else
						IND=org.t_strial(i):org.t_etrial(i);
					end
					newtheta(IND)=theta(IND);
					costheta(IND)=cos(theta(IND));
					sintheta(IND)=sin(theta(IND));
				end
				org.sp_t(:,1)=newtheta;
				% org.cossp_t(:,1)=costheta;
				% org.sinsp_t(:,1)=sintheta;
			end
			function org=fix_Rot_motion_org(org)
				for i=1:numel(org.t_strial)
					if org.motionLog(i)==1
						if i~=numel(org.t_strial)
							theta=org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1);
						else
							theta=org.sp_t(org.t_strial(i):org.t_etrial(i),1);
						end
						theta=force_monotonic_increase(theta,min=deg2rad(45/2),max=deg2rad(360+45/2));				
						if i~=numel(org.t_strial)
							org.sp_t(org.t_strial(i):org.t_strial(i+1)-1,1)=theta;
						else
							org.sp_t(org.t_strial(i):org.t_etrial(i),1)=theta;
						end
					end
				end
			end
			function org=fixROT(org)
				theta=org.sp_t(:,1);
				ind=find(org.mtypes=="rot");
				for i=1:length(org.t_strial)
					if org.motionLog(i)~=ind
						continue;
					end
					if i<length(org.t_strial)
						IND=org.t_strial(i):org.t_strial(i+1);
					else
						IND=org.t_strial(i):org.t_etrial(i);
					end
					tpk=theta(IND);
					tpk=atan2(org.sinsp_t(IND,1),org.cossp_t(IND,1));
					% 
					% % Assume tpk is [T, P, N] (time x parameters x trials)
					% T = size(tpk,1);
					% N = size(tpk,3);
					% 
					% % 1. Determine where all parameters are increasing between consecutive time points.
					% inc = squeeze(all(diff(tpk, 1, 1) > 0, 2));  % Now inc is [T-1, N]
					% 
					% % 2. Append false at the beginning and end for each trial, then compute diff to mark patch boundaries.
					% inc_ext = [false(1,N); inc; false(1,N)];  % now size: [T+1, N]
					% d = diff(double(inc_ext));                % d==1 marks start, d==-1 marks end; size: [T, N]
					% 
					% % 3. For each trial, find start and end indices of all increasing patches.
					% startCells = arrayfun(@(j) find(d(:,j)==1), 1:N, 'UniformOutput', false);
					% endCells   = arrayfun(@(j) find(d(:,j)==-1)-1, 1:N, 'UniformOutput', false);
					% 
					% % For each trial, choose the longest patch.
					% [~,idxLongest] = cellfun(@(s,e) max(e-s+1), startCells, endCells);
					% startTime = cellfun(@(s,i) s(i), startCells, num2cell(idxLongest));
					% endTime   = cellfun(@(e,i) e(i)+1, endCells, num2cell(idxLongest)); % add one because diff is one less than tpk
					% 
					% % 4. Create a logical mask that is true only for the longest patch in each trial.
					% tvec = (1:T)';  % column vector of time indices
					% mask = (tvec >= startTime) & (tvec <= endTime);  % uses implicit expansion; size: [T, N]
					% mask = reshape(mask, T, 1, N);                   % reshape to 3D: [T, 1, N]
					% mask = repmat(mask, [1, size(tpk,2), 1]);         % replicate across the parameter dimension
					% 
					% % 5. Set values outside the longest monotonic increasing patch to NaN.
					% tpk(~mask) = NaN;
					% tpk=fillmissing(tpk,"previous",1);
					% tpk=fillmissing(tpk,"next",1);
					% 
					% if any(isnan(tpk),"all"); error("There are NaN values in tpk");end
					% if any(diff(tpk,1,1)<0)
					% 	warning("There still remains timepoints where phi is not monotonically increasing");
					% end
					theta(IND)=tpk;
				end
				org.sp_t(:,1)=theta;
			end
			function sp_mtype=fixRotMotionSP(sp_mtype,org)
				ind=org.ind_mtypes("rot");
				tpk=sp_mtype{ind}(:,1,:);

				% Assume tpk is [T, P, N] (time x parameters x trials)
				T = size(tpk,1);
				N = size(tpk,3);

				% 1. Determine where all parameters are increasing between consecutive time points.
				inc = squeeze(all(diff(tpk, 1, 1) > 0, 2));  % Now inc is [T-1, N]

				% 2. Append false at the beginning and end for each trial, then compute diff to mark patch boundaries.
				inc_ext = [false(1,N); inc; false(1,N)];  % now size: [T+1, N]
				d = diff(double(inc_ext));                % d==1 marks start, d==-1 marks end; size: [T, N]

				% 3. For each trial, find start and end indices of all increasing patches.
				startCells = arrayfun(@(j) find(d(:,j)==1), 1:N, 'UniformOutput', false);
				endCells   = arrayfun(@(j) find(d(:,j)==-1)-1, 1:N, 'UniformOutput', false);

				% For each trial, choose the longest patch.
				[~, idxLongest] = cellfun(@(s,e) max(e-s+1), startCells, endCells);
				startTime = cellfun(@(s,i) s(i), startCells, num2cell(idxLongest));
				endTime   = cellfun(@(e,i) e(i)+1, endCells, num2cell(idxLongest)); % add one because diff is one less than tpk

				% 4. Create a logical mask that is true only for the longest patch in each trial.
				tvec = (1:T)';  % column vector of time indices
				mask = (tvec >= startTime) & (tvec <= endTime);  % uses implicit expansion; size: [T, N]
				mask = reshape(mask, T, 1, N);                   % reshape to 3D: [T, 1, N]
				mask = repmat(mask, [1, size(tpk,2), 1]);         % replicate across the parameter dimension

				% 5. Set values outside the longest monotonic increasing patch to NaN.
				tpk(~mask) = NaN;
				tpk=fillmissing(tpk,"previous",1);
				tpk=fillmissing(tpk,"next",1);

				if any(isnan(tpk),"all"); error("There are NaN values in tpk");end
				if any(diff(tpk,1,1)<0)
					warning("There still remains timepoints where phi is not monotonically increasing");
				end
				sp_mtype{ind}(:,1,:)=tpk;
			end
			function cellac_downsampled=Downsample_cellac(cellac,timestamp,org)
				x=1:size(cellac,1);
				xq=org.f2p_t(timestamp);
				cellac_downsampled=interp1(x,cellac,xq);
				if any(isnan(cellac_downsampled),"all")
					error("NaN values detected in cellac_downsampled")
				end
			end
			function DLC_downsampled=Downsample_DLC(DLC_raw,timestamp,org)
				fdlc_ind=org.fdlc_t(timestamp);
				fdlc_ind=round(fdlc_ind);
				DLC_downsampled=applyFunctionLayers(DLC_raw,@(a)a(fdlc_ind,:));
			end
			function trials=Downsample_trial_border(org,timestamp)
				trials=[];
				trials.t_strial=interp1(timestamp,[1:numel(timestamp)]',org.t_strial,"nearest");
				trials.t_etrial=interp1(timestamp,[1:numel(timestamp)]',org.t_etrial,"nearest");
				trials.mtypes=org.mtypes;
				trials.names_trial=strings(size(org.motionLog));
				for i=1:length(org.mtypes)
					ind=org.motionLog==i;
					trials.names_trial(ind)=org.mtypes(i);
				end
			end
			function DLC=DLC_interpolateHighErrorDataPoints(DLC,triErrorThresh)
				Qind1=cellfun(@(a)size(a,1),DLC.tXYZbp);
				Qind2=cellfun(@(a)size(a,2),DLC.tXYZbp);
				Q=cell2mat(DLC.tXYZbp(:));
				triError=DLC.triError(:);
				lowErrorlog=triError<=triErrorThresh;
				% make sure that the endpoints of each body part do not get
				% mixed up in interpolation
				for i=1:length(DLC.tXYZbp)
					ind=(i-1)*Qind1(1)+1:i*Qind1(1);
					a=Q(ind,:);
					lowE=lowErrorlog(ind);
					if lowE(1)==0
						a(1,:)=a(find(lowE,1,"first"),:);
						lowE(1)=1;
					end
					if lowE(end)==0
						a(end,:)=a(find(lowE,1,"last"),:);
						lowE(end)=1;
					end
					lowErrorlog(ind)=lowE;
					Q(ind,:)=a;
				end
				x=find(lowErrorlog);
				y=Q(lowErrorlog,:);
				xq=[1:size(Q,1)]';
				interpy=interp1(x,y,xq);
				interpy=fillmissing(interpy,"previous",1);
				interpy=fillmissing(interpy,"next",1);
				newQ=mat2cell(interpy,Qind1,size(interpy,2));
				newQ=reshape(newQ,size(DLC.tXYZbp));
				DLC.tXYZbp=newQ;
			end
			function [DLCang,DLCdot]=DLC_make_Joint_Angle(Q_bp,rec3d)
				DLCang=[];
				DLCdot=[];
				for i=1:size(rec3d.jointList,1)
					inds=rec3d.ind_bp(rec3d.jointList(i,:));
					Q=Q_bp(inds);
					a=Q{2}-Q{1};	a=a./sqrt(dot(a,a,2));
					b=Q{3}-Q{1};	b=b./sqrt(dot(b,b,2));
					ang=acos(dot(a,b,2));
					DLCang(:,i)=ang;
					DLCdot(:,i)=dot(a,b,2);
				end
			end
			function DLC=DLC_calculateDerivatives(DLC,period_ms_cutoff,DivideDerivWithTime,org)
				DLC.tDXYZbp=[];
				DLC.tD2XYZbp=[];
				DLCsamplingRate=1000*numel(org.t_fdlc)./(org.t_fdlc(end)-org.t_fdlc(1));
				for i=1:length(DLC.tXYZbp)
					[DLC.tDXYZbp{i},DLC.tD2XYZbp{i},filtered_data]=calculateDerivatives_along1stdim(DLC.tXYZbp{i},period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate);
					%%%%%%%% plot filtered data
					if false
						figure;plot(DLC.tXYZbp{i}(:,1),LineWidth=2);hold on;plot(filtered_data(:,1),LineWidth=2);
						xticks(0:1000:size(filtered_data,1))
						xlabel("time (ms)")
						legend("original data","filtered data")
						title("all the periods shorter than 100 ms is removed")
						ChangeFigureSettings(TitleSize=20,XLabelSize=15,LegendSize=15)
					end
				end
				[DLC.tDANGbp,DLC.tD2ANGbp,filtered_data]=calculateDerivatives_along1stdim(DLC.tANGbp,period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate);
				[DLC.tDDOTbp,DLC.tD2DOTbp,filtered_data]=calculateDerivatives_along1stdim(DLC.tDOTbp,period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate);
			end
			function smoothed_DLC=smooth_DLC(DLC,period_cutoff_sec)
				if ~isempty(period_cutoff_sec)
					smoothed_DLC=DLC;
					for i=1:length(DLC.tXYZbp)
						smoothed_DLC.tXYZbp{i}=lowpass_filter(DLC.tXYZbp{i},1000,1/period_cutoff_sec);
					end
					smoothed_DLC.tANGbp=lowpass_filter(DLC.tANGbp,1000,1/period_cutoff_sec);
					smoothed_DLC.tDOTbp=lowpass_filter(DLC.tDOTbp,1000,1/period_cutoff_sec);
				else
					smoothed_DLC=DLC;
				end
			end
			function [ddata,d2data,filtered_data]=calculateDerivatives_along1stdim(data,period_ms_cutoff,DivideDerivWithTime,DLCsamplingRate)
				cutoff=floor(size(data,1)/period_ms_cutoff)+1;
				% filtered_data=lowpass_filter_mat_along1stdim(data,cutoff,Display_WarningIfNonReal=false);
				filtered_data=lowpass_filter(data,300,1000/period_ms_cutoff);	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sampling rate is 300, not 1000	2025/06/29
				ddata=diff(filtered_data,1,1);
				if DivideDerivWithTime
					ddata=ddata./(1/DLCsamplingRate);
				end

				% changed acceleration derivation 2025/06/29 %%%%%%%%%%%%%%
				% d2data=diff(filtered_data,2,1);
				sampling_rate=DLCsamplingRate;
				diff_sec=0.1;
				diff_time_points=round(sampling_rate*diff_sec);
				d2data=ddata(diff_time_points:end,:)-ddata(1:end-diff_time_points+1,:);
				d2data(size(d2data,1)+1:size(ddata,1),:)=NaN;
				d2data=fillmissing(d2data,"previous",1);
				if DivideDerivWithTime
					d2data=d2data./(diff_sec);
				end
			end
		end
		%%
		function DLC_raw=make_modified_DLC_raw(obj,DLC_numeric,rec3d)
			Q_bp=mat2cell(DLC_numeric,size(DLC_numeric,1),repmat(3,1,size(DLC_numeric,2)/3));
			DLC=[];
			DLC.tXYZbp					=Q_bp;
			DLC.triError				=cell2mat(rec3d.offset(:)');
			DLC=DLC_interpolateHighErrorDataPoints(DLC,obj.triErrorThresh);
			[DLC.tANGbp,DLC.tDOTbp]		=DLC_make_Joint_Angle(Q_bp,rec3d);
			period_ms_cutoff=round(obj.PeriodCutoff4Derivatives_sec*1000);
			DLC							=DLC_calculateDerivatives(DLC,period_ms_cutoff);
			DLC							=smooth_DLC(DLC,obj.PeriodCutoff4Derivatives_sec);
			DLC_raw=DLC;
			function DLC=DLC_interpolateHighErrorDataPoints(DLC,triErrorThresh)
				Qind1=cellfun(@(a)size(a,1),DLC.tXYZbp);
				Qind2=cellfun(@(a)size(a,2),DLC.tXYZbp);
				Q=cell2mat(DLC.tXYZbp(:));
				triError=DLC.triError(:);
				lowErrorlog=triError<=triErrorThresh;
				% make sure that the endpoints of each body part do not get
				% mixed up in interpolation
				for i=1:length(DLC.tXYZbp)
					ind=(i-1)*Qind1(1)+1:i*Qind1(1);
					a=Q(ind,:);
					lowE=lowErrorlog(ind);
					if lowE(1)==0
						a(1,:)=a(find(lowE,1,"first"),:);
						lowE(1)=1;
					end
					if lowE(end)==0
						a(end,:)=a(find(lowE,1,"last"),:);
						lowE(end)=1;
					end
					lowErrorlog(ind)=lowE;
					Q(ind,:)=a;
				end
				x=find(lowErrorlog);
				y=Q(lowErrorlog,:);
				xq=[1:size(Q,1)]';
				interpy=interp1(x,y,xq);
				interpy=fillmissing(interpy,"previous",1);
				interpy=fillmissing(interpy,"next",1);
				newQ=mat2cell(interpy,Qind1,size(interpy,2));
				newQ=reshape(newQ,size(DLC.tXYZbp));
				DLC.tXYZbp=newQ;
			end
			function [DLCang,DLCdot]=DLC_make_Joint_Angle(Q_bp,rec3d)
				DLCang=[];
				DLCdot=[];
				for i=1:size(rec3d.jointList,1)
					inds=rec3d.ind_bp(rec3d.jointList(i,:));
					Q=Q_bp(inds);
					a=Q{2}-Q{1};	a=a./sqrt(dot(a,a,2));
					b=Q{3}-Q{1};	b=b./sqrt(dot(b,b,2));
					ang=acos(dot(a,b,2));
					DLCang(:,i)=ang;
					DLCdot(:,i)=dot(a,b,2);
				end
			end
			function DLC=DLC_calculateDerivatives(DLC,period_ms_cutoff)
				DLC.tDXYZbp=[];
				DLC.tD2XYZbp=[];
				for i=1:length(DLC.tXYZbp)
					[DLC.tDXYZbp{i},DLC.tD2XYZbp{i},filtered_data]=calculateDerivatives_along1stdim(DLC.tXYZbp{i},period_ms_cutoff);
					%%%%%%%% plot filtered data
					if false
						figure;plot(DLC.tXYZbp{i}(:,1),LineWidth=2);hold on;plot(filtered_data(:,1),LineWidth=2);
						xticks(0:1000:size(filtered_data,1))
						xlabel("time (ms)")
						legend("original data","filtered data")
						title("all the periods shorter than 100 ms is removed")
						ChangeFigureSettings(TitleSize=20,XLabelSize=15,LegendSize=15)
					end
				end
				[DLC.tDANGbp,DLC.tD2ANGbp,filtered_data]=calculateDerivatives_along1stdim(DLC.tANGbp,period_ms_cutoff);
				[DLC.tDDOTbp,DLC.tD2DOTbp,filtered_data]=calculateDerivatives_along1stdim(DLC.tDOTbp,period_ms_cutoff);
			end
			function smoothed_DLC=smooth_DLC(DLC,period_cutoff_sec)
				if ~isempty(period_cutoff_sec)
					smoothed_DLC=DLC;
					for i=1:length(DLC.tXYZbp)
						smoothed_DLC.tXYZbp{i}=lowpass_filter(DLC.tXYZbp{i},1000,1/period_cutoff_sec);
					end
					smoothed_DLC.tANGbp=lowpass_filter(DLC.tANGbp,1000,1/period_cutoff_sec);
					smoothed_DLC.tDOTbp=lowpass_filter(DLC.tDOTbp,1000,1/period_cutoff_sec);
				else
					smoothed_DLC=DLC;
				end
			end
			function [ddata,d2data,filtered_data]=calculateDerivatives_along1stdim(data,period_ms_cutoff)
				cutoff=floor(size(data,1)/period_ms_cutoff)+1;
				% filtered_data=lowpass_filter_mat_along1stdim(data,cutoff,Display_WarningIfNonReal=false);
				filtered_data=lowpass_filter(data,1000,1000/period_ms_cutoff);
				ddata=diff(filtered_data,1,1);
				d2data=diff(filtered_data,2,1);
			end
		end
		%%
		function setCurrentData_to_continuous(obj)
			arguments
				obj
			end
			obj.tck = obj.cellac_continuous;
			obj.tpk = obj.sp_continuous;
			obj.DLC = obj.DLC_countinuous;
			obj.invDLC=obj.invDLC_countinuous;


			obj.make_DLC_speed();
			obj.make_DLC_accel();

			% Preprocess the SP data.
			obj.tpk = PreprocessSPdata(obj.tpk);

			% obtain T
			obj.SP_T=SPparams2TransMat(obj.tpk);

			% calculate DLC inverse transformed by stewart platform
			% obj.invDLC = applyInverseSPTrans2DLC(obj.DLC, obj.SP_T,1/obj.SamplingRate_Hz);

			% obtian all numeric DLC data (still a struct though)
			obj.DLC_allNumeric=obj.DLC;
			DLCflds=string(fields(obj.DLC_allNumeric));
			for k=DLCflds(:)'
				try
					obj.DLC_allNumeric.(k)=cell2mat(obj.DLC_allNumeric.(k));
				end
			end

			% obtain all numeric invDLC data (still a struct though)
			obj.invDLC_allNumeric=obj.invDLC;
			invDLCflds=string(fields(obj.invDLC_allNumeric));
			for k=invDLCflds(:)'
				try
					obj.invDLC_allNumeric.(k)=cell2mat(obj.invDLC_allNumeric.(k));
				end
			end	% concatenate data

			obj.DLC_inv_or_normal="normalPosture";
		end
		%% loadData_with_trial_margin
		function loadData_with_trial_margin(obj,options)
			arguments
				obj
				options.pre_margin_in_sec=0
				options.post_margin_in_sec=0
			end
			if isempty(obj.org)
				obj.loadData();
			end
			premrgn=options.pre_margin_in_sec*1000;
			postmrgn=options.post_margin_in_sec*1000;
			org=obj.org;
			prevorg=org;
			org.t_strial=org.t_strial-premrgn;
			org.t_etrial=org.t_etrial+postmrgn;
			validtrials=find(org.t_strial>=prevorg.t_strial(1) & org.t_etrial<=prevorg.t_etrial(end));
			org.t_strial=org.t_strial(validtrials);
			org.t_etrial=org.t_etrial(validtrials);
			org.motionLog=org.motionLog(validtrials);
			org.strial_t=[];
			org.etrial_t=[];
			obj.loadData(UseCachedData=true,Cached_org=org);
		end
	end
	properties
		full_rec3d
		DLCvid
	end
	methods
		function loadFullDLC(obj)
			Dirlist=load("/path/to/data/AllMiceData/ListOf_original_DirectoriesForAllMice.mat");
			allmicedir=Dirlist.allmicedir;
			temp=load(allmicedir{obj.plane}.DLCdata);
			obj.full_rec3d=temp.rec3d;
		end
	end
	properties (Constant)
		LimbNames_for_DLCvid= ["RightHand1", "RightElbow1", "RightToe1", "RightHeel1", "RightKnee1", ...
            "BackU1", "BackM1", "BackL1", "Tail1", ...
            "LeftHand2", "LeftElbow2", "LeftToe2", "LeftHeel2", "LeftKnee2", ...
            "BackU2", "BackM2", "BackL2", "Tail2", ...
            "LeftHand3", "LeftElbow3", "LeftToe3", "LeftHeel3", "LeftKnee3", ...
            "RightHand4", "RightElbow4", "RightToe4", "RightHeel4", "RightKnee4"]
	end
	methods
		%% DLC video, execute after continuous data load
		function load_DLCVideo(obj)
			obj.DLCvid=VideoReader(obj.rec3d.DLCvid);
		end
		function fdlc_ind=get_DLCvid_frameIND(obj)
			timestamp=obj.timestamp;
			org=obj.org;
			fdlc_ind=org.fdlc_t(timestamp);
			fdlc_ind=round(fdlc_ind);
		end
		function Image=get_DLCFrame(obj,DLCtimepoints)
			fdlc_ind=obj.get_DLCvid_frameIND();
			Image=read(obj.DLCvid,fdlc_ind(DLCtimepoints));
		end
		
		function XYbp=get_2d_bodyPart_for_DLC(obj)
			XYbp=obj.full_rec3d.imageCoord_bp;
		end
		function indicesCell=get_DLCind_for_body_parts(obj)
			Limbs = obj.LimbNames_for_DLCvid;
			limbList = convertStringsToChars(Limbs);
			cameraNumbers = cellfun(@(x) str2double(regexp(x, '\d+$', 'match', 'once')), limbList);
			
			% Remove trailing digits to obtain base limb names.
			baseNames = cellfun(@(x) regexprep(x, '\d+$', ''), limbList, 'UniformOutput', false);

			% Get unique base names in the order of first appearance.
			[uniqueBaseNames, ~, ~] = unique(baseNames, 'stable');
			uniqueBaseNames=convertCharsToStrings(uniqueBaseNames);

			% Convert the unique names to a string array.
			limbNames = string(uniqueBaseNames);

			% Create a cell array where each element contains the indices of limbList
			% corresponding to each unique limb name.
			indicesCell = cell(size(uniqueBaseNames));
			for i = 1:length(uniqueBaseNames)
				indicesCell{i} = find(strcmp(baseNames, uniqueBaseNames(i)));
			end
		end
		function cam_tXY_bpbp=get_camera_numbers(obj)
			arguments
				obj
			end
			Limbs = obj.LimbNames_for_DLCvid;
			limbList = convertStringsToChars(Limbs);
			cameraNumbers = cellfun(@(x) str2double(regexp(x, '\d+$', 'match', 'once')), limbList);
			
			% Remove trailing digits to obtain base limb names.
			baseNames = cellfun(@(x) regexprep(x, '\d+$', ''), limbList, 'UniformOutput', false);

			% Get unique base names in the order of first appearance.
			[uniqueBaseNames, ~, ~] = unique(baseNames, 'stable');
			uniqueBaseNames=convertCharsToStrings(uniqueBaseNames);

			% Convert the unique names to a string array.
			limbNames = string(uniqueBaseNames);

			% Create a cell array where each element contains the indices of limbList
			% corresponding to each unique limb name.
			indicesCell = cell(size(uniqueBaseNames));
			for i = 1:length(uniqueBaseNames)
				indicesCell{i} = find(strcmp(baseNames, uniqueBaseNames(i)));
			end

			% Extract the camera number from the end of each string in limbList.
			cameraNumbers = cellfun(@(x) str2double(regexp(x, '\d+$', 'match', 'once')), limbList);

			dlc_tXY_bpbp=[];
			cam_tXY_bpbp=[];
			for i=1:length(indicesCell)
				cam_tXY_bpbp{i}=cameraNumbers(indicesCell{i});
			end
		end
		function overlayBodyEdgesOnDLCimage(obj,timepoint)
			arguments
				obj
				timepoint
			end
			currframe=obj.get_DLCvid_frameIND();
			currframe=currframe(timepoint);
			bp3d=obj.DLC_countinuous.tXYZbp;
			bp2d=obj.get_2d_bodyPart_for_DLC();
			Numof2Ddata=sum(cellfun(@(a)numel(a),bp2d));
			% plot DLC image
			catbp2d=cat(2,bp2d{:});
			catbp2d=cellfun(@(a)a(currframe,:),catbp2d,UniformOutput=false);
			catbp2d=cat(1,catbp2d{:});
			cmap_bp2d=make_color_4_bp2d();
			scatter(catbp2d(:,1),catbp2d(:,2),10,cmap_bp2d,"filled");

			% get bp ind
			bp2d_ind=obj.get_DLCind_for_body_parts();

			% get indices for each bodypart
			bp2d_bpind=[];
			for i=1:numel(obj.rec3d.bp_ind)
				bp=obj.rec3d.bp_ind(i);
				ind=find(containsSubstring(obj.LimbNames_for_DLCvid,bp));
				bp2d_bpind(i)=find(cellfun(@(a)all(ismember(ind,a),"all"),bp2d_ind));
			end

			camind_bp2d=get_camera_numbers();

			% connect bp
			for i=1:size(obj.rec3d,edgeList)
				bps=csp.rec3d.edgeList(i,:);
				bpind=obj.rec3d.ind_bp(bps);
				cam=[];
				cam{1}=camind_bp2d{bp2d_bpind(bpind(1))};
				cam{2}=camind_bp2d{bp2d_bpind(bpind(2))};
				for c1=cam{1}(:)'
					for c2=cam{2}(:)'
						
					end
				end
			end
			function cmap=make_color_4_bp2d()
				colors=hsv(numel(bp2d));
				cmap=repelem(colors,2,1);
			end
		end
	end
	%% Atlas Image modification
	methods
		function inv_image=inverse_AtlasImage(obj)
			thresh=50;
			inv_image=255-obj.AtlasImage;
			inv_image(inv_image<=thresh)=0;
			inv_image(inv_image>=thresh)=255;
			inv_image(:,:,[1,3])=0;	% set to green
		end
		function image=normal_AtlasImage(obj)
			inv_image=obj.inverse_AtlasImage();
			image=255-inv_image;
		end
		function mask=make_AtlasImage_mask(obj)
			im=obj.AtlasImage;
			Log=prod(im>=200,3);
			mask=~Log;
		end
		function inv_image=inverse_HalfAtlas(obj)
			thresh=50;
			inv_image=255-obj.HalfAtlas;
			inv_image(inv_image<=thresh)=0;
			inv_image(inv_image>=thresh)=255;
			inv_image(:,:,[1,3])=0;	% set to green
		end
		function image=normal_HalfAtlas(obj)
			inv_image=obj.inverse_HalfAtlas();
			image=255-inv_image;
		end
		function mask=make_HalfAtlas_mask(obj)
			im=obj.HalfAtlas;
			Log=prod(im>=200,3);
			mask=~Log;
		end
		%%
		function makeDirectoryList(obj,options)
			arguments
				obj
				options.coordinate (1,:){mustBeNumeric}=[]
				options.filename=[]
				options.dimNames=[]
				options.label=[]
			end
			IND=num2cell(options.coordinate);
			if ~isfield(obj.savedDirectories,options.label)
				obj.savedDirectories.(options.label)=string;
			end
			obj.savedDirectories.(options.label)(IND{:})=options.filename;
			if ~isempty(options.dimNames)
				obj.savedDirectories.dimNames=options.dimNames;
			end
		end
		% function clearSaveAndMakeDirectoryList(obj,options)
		% 	arguments
			% 	obj
			% 	options.Properties2Erase (1,:) string = []
			% 	options.Properties2keep (1,:) string = []
			% 	options.filename =[]
		% 
			% 	options.coordinate (1,:){mustBeNumeric}=[]
			% 	options.dimNames=[]
			% 	options.label=[]
		% 	end
		% 	options1 = struct( ...
			% 	"Properties2Erase", options.Properties2Erase, ...
			% 	"Properties2Keep",  options.Properties2Keep, ...
			% 	"filename",         options.filename ...
			% 	);
		% 	options1=namedargs2cell(options1);
		% 	options2 = struct( ...
			% 	"coordinate", options.coordinate, ...
			% 	"filename",  options.filename, ...
			% 	"dimNames",   options.dimNames, ...
			% 	"label",      options.label ...
			% 	);
		% 	options2=namedargs2cell(options2);
		% 	obj.clearAndSave(options1{:});
		% 	obj.makeDirectoryList(options2{:});
		% end
		%%
	end
	%% make a flipped image
	properties
		x_mean
		HalfAtlas
		HalfROI
	end
	methods
		function calculatemidline(obj,options)
			arguments
				obj
				options.disp_data (1,1) logical =false
			end
			leftreg=obj.contour_LeftRegion;
			rightreg=obj.contour_RightRegion;
			leftmean=getMeanCoord_AllContours(leftreg);
			rightmean=getMeanCoord_AllContours(rightreg);
			obj.x_mean=mean((leftmean(:,2)+rightmean(:,2))/2);
			if options.disp_data
				disp("y diff")
				leftmean(:,1)-rightmean(:,1)
				disp("x mean")
				(leftmean(:,2)+rightmean(:,2))/2
			end
			function meanCoord=getMeanCoord_AllContours(cellarray)
				meanCoord=[];
				for i=1:numel(cellarray)
					meanCoord(i,:)=getMeanCoord(cellarray{i});
				end
			end
			function meanCoord=getMeanCoord(logmat)
				% a is your logical 2D array
				[rows, cols] = find(logmat);        % get row‑ and column‑indices of all true pixels
				meanRow    = mean(rows);       % average row index
				meanCol    = mean(cols);       % average column index
				meanCoord  = [meanRow, meanCol];  % [y, x] format
			end
		end
		function makeHalfAtlas(obj)
			arguments
				obj
			end
			obj.calculatemidline();
			obj.HalfAtlas=obj.AtlasImage(:,round(obj.x_mean):end,:);
		end
		function ROI=makeHalfCoord(obj,ROI)
			arguments
				obj
				ROI=[]
			end
			emptyROI=false;
			if isempty(ROI)
				emptyROI=true;
				ROI=obj.regROI;
			end
			if size(ROI)~=2
				error("The size in dim 2 of ROI has to be 2")
			end
			% regROI(:,2) is the x coordinate and regROI(:,1) is the y coordinate
			% (Vertically down) the origin is top left of image
			obj.calculatemidline();
			obj.makeHalfAtlas();
			x=ROI(:,2);
			x=x-obj.x_mean+1;
			negxlog=x<0;
			x(negxlog)=abs(x(negxlog))+1;
			ROI(:,2)=x;
			if emptyROI
				obj.HalfROI=ROI;
			end
		end
	end
	%% storeProperties2StructCell
	properties
		storedProperties
		storedProperties_DimValueList
	end
	methods
		function storeProperties2StructCell(obj,options)
			%{
    		Function: storeProperties2StructCell
    		--------------------------------------
    		This function stores selected properties from the input object 'obj' into a structured cell array
    		based on dimension-specific options provided in the 'options' structure. The options include:
      		- 'Properties2Store': a list of property names from 'obj' to be stored.
      		- 'Properties4StructCellDim': a list of property names that define dimensions for the structured cell.
      		- 'ValueListOf_Properties4StructCellDim': reference lists for each dimension property.
    		
    		The process involves:
      		1. Obtaining an initial structured cell (or reusing an existing one).
      		2. Determining indices for each dimension using the given reference values.
      		3. Assigning the corresponding properties from 'obj' into the structured cell at the computed indices.
      		4. Optionally storing the dimension value lists separately in 'obj.storedProperties_DimValueList'
         		if the dimension properties and corresponding values are provided.
    		
    		Helper nested functions:
      		- obtainInitialStructCell: Retrieves or initializes the structured cell array.
      		- getIND4EachProperties: Computes the index for each of the dimension properties.
      		- InitializeStruct: (Unused in the current implementation) Initializes a structure array with empty fields.
    		
    		Example Usage:
      		csp.storeProperties2StructCell(...
				Properties2Store=["Meanscore_CellRidgeparams","Stdscore_CellRidgeparams","ridgeparams"], ...
           		Properties4StructCellDim=["nameOfCurrentSettingRidge","TimeShiftInSec"], ...
           		ValueListOf_Properties4StructCellDim={["SP_and_DLC","SP_only","DLC_only"], shiftInSec})
			%}
			arguments
				obj
				options.Properties2Store			(1,:) string
				options.Properties4StructCellDim	(1,:) string
				options.ValueListOf_Properties4StructCellDim	(1,:) cell
			end
			data=obtainInitialStructCell(options);
			INDcell=getIND4EachProperties(options.Properties4StructCellDim,options.ValueListOf_Properties4StructCellDim);
			for prop=options.Properties2Store(:)'
				data.(prop){INDcell{:}}=obj.(prop);
			end
			obj.storedProperties=data;
			if ~isempty(options.Properties4StructCellDim) && ~isempty(options.ValueListOf_Properties4StructCellDim)
				obj.storedProperties_DimValueList=[];
				for i=1:length(options.Properties4StructCellDim)
					prop=options.Properties4StructCellDim(i);
					obj.storedProperties_DimValueList.(prop)=options.ValueListOf_Properties4StructCellDim{i};
				end
			end
			function INDcell=getIND4EachProperties(PropertyNames,ReferenceValues4EachProp)
				INDcell={};
				for k=1:numel(PropertyNames)
					INDcell{1,k}=find(ReferenceValues4EachProp{k}==obj.(PropertyNames(k)));
				end
			end
			function data=obtainInitialStructCell(options)
				if isempty(obj.storedProperties)
					sz=cellfun(@(a)numel(a),options.ValueListOf_Properties4StructCellDim);
					data=struct;
					for prop=options.Properties2Store(:)'
						if numel(sz)==1
							data.(prop)=cell([sz,1]);
						else
							data.(prop)=cell(sz);
						end
					end
				else
					data=obj.storedProperties;
				end
			end
			% function data = InitializeStruct(SizeOfStructArray, ListOfFields)
			% 	data = struct;
			% 	ind = num2cell(SizeOfStructArray);
			% 	for k = ListOfFields(:)'
				% 	data(ind{:}).(k) = [];
			% 	end
			% end
		end
		%% example use
		% csp=CSP();
		% csp.loadData(1);
		% shiftInSec=[-1:0.2:1];
		% csp.TimeShifted_Initialize_timeshift(TimeShiftInSec_vec=shiftInSec);
		% for i=1:length(shiftInSec)
		% 	csp.makeTimeShiftedData(TimeShiftInSec=shiftInSec(i));
		% 	csp.storeProperties2StructCell(Properties2Store=["p8_tckd"], ...
			% 	Properties4StructCellDim=["currPlaneName4Directory","TimeShiftInSec"], ...
			% 	ValueListOf_Properties4StructCellDim= ...
			% 	{csp.PlaneName4Directory,shiftInSec});
		% 	i
		% end
		% csp.figureTile;
		% for i=1:length(shiftInSec)
		% 	nexttile;
		% 	plot(csp.storedProperties.p8_tckd{1,i}(:,1,1,1));
		% end
		% ChangeFigureSettings(Linkaxes="xy")
	end
	%%
	properties
		FilePathsList
		FilePathsList_Labels4EachCategory
	end
	methods
		function OrganizePathsIntoCategories(obj,options)
			arguments
				obj
				options.FolderName (1,1) string
				options.strSeparator (1,:) string
				options.Categories (1,:) string
				options.NumericCategoryIndex (1,:)=[]
				options.FinalCellDimFormat (1,:)
				options.AllLabels4EachCategory cell={}
			end
			if length(options.strSeparator)+1~=length(options.Categories);error("Length does not match");end
			categories=options.Categories;
			valid_categories=categories(categories~="");
			Allfiles=getAllFiles(options.FolderName);
			listcell=splitString2Chunks(Allfiles,options.strSeparator,direction="last");
			[LabelIND_file_4eachCategory,LabelVal_4eachCategory]=makeList4EachCategory(listcell,categories,options.NumericCategoryIndex);
			[LabelIND_file_4eachCategory,ValidFileLog]=removeNaNfiles_fromAllCategories(LabelIND_file_4eachCategory);
			cellIND_4eachFile=storeArraysInStruct2SpecifiedOrder(LabelIND_file_4eachCategory,options.FinalCellDimFormat);
			sizsOfFinalCellArray=obtainLengthInStruct2SpecifiedOrder(LabelVal_4eachCategory,options.FinalCellDimFormat);
			IND=sub2ind(sizsOfFinalCellArray,cellIND_4eachFile{:});
			AllDirList=strings(sizsOfFinalCellArray);
			AllDirList(IND)=Allfiles(ValidFileLog);
			obj.FilePathsList=AllDirList;
			obj.FilePathsList_Labels4EachCategory=LabelVal_4eachCategory;

			function sz=obtainLengthInStruct2SpecifiedOrder(Struct,OrderStr)
				sz=[];
				for i=OrderStr(:)'
					sz(end+1)=length(Struct.(i));
				end
			end
			function celldata=storeArraysInStruct2SpecifiedOrder(Struct,OrderStr)
				celldata={};
				for i=OrderStr(:)'
					celldata{end+1}=Struct.(i)(:);
				end
			end
			function [Struct,NoNanLog]=removeNaNfiles_fromAllCategories(Struct)
				flds=string(fields(Struct));
				catdata=[];
				for i=1:length(flds)
					catdata=cat(2,catdata,Struct.(flds(i))(:));
				end
				NoNanLog=logical(prod(~isnan(catdata),2));
				Struct=structfun(@(a)a(NoNanLog),Struct,UniformOutput=false);
			end
			function [LabelIND_file,LabelVal]=makeList4EachCategory(listcell,categories,NumericCategoryIndex)
				LabelIND_file=struct;
				LabelVal=struct;
				for i=1:length(categories)
					if categories(i)==""
						continue;
					end
					if ismember(i,NumericCategoryIndex)
						list4eachFile=cellfun(@(strArray)IndexStrArray_convert2Num(strArray,i),listcell,UniformOutput=true);
						if isempty(options.AllLabels4EachCategory)
							uniqueLabels=unique(list4eachFile);
							uniqueLabels=uniqueLabels(~isnan(uniqueLabels));
							uniqueLabels=sort(uniqueLabels);
						else
							uniqueLabels=options.AllLabels4EachCategory{i};
						end
					else
						list4eachFile=cellfun(@(strAr)IndexStrArray(strAr,i),listcell,UniformOutput=true);
						if isempty(options.AllLabels4EachCategory)
							uniqueLabels=unique(list4eachFile);
							uniqueLabels=uniqueLabels(uniqueLabels~="");
						else
							uniqueLabels=options.AllLabels4EachCategory{i};
						end
					end
					LabelVal.(categories(i))=uniqueLabels;
					LabelIND_file.(categories(i))=arrayfun(@(a)find_but_NaN_if_empty(uniqueLabels==a),list4eachFile);
					checkLength(list4eachFile,listcell);
				end
				function ind=find_but_NaN_if_empty(logicalArray)
					ind=find(logicalArray);
					if isempty(ind)
						ind=NaN;
					end
				end
				function str=IndexStrArray(stringArray,ind)
					str="";
					try
						str=stringArray(ind);	  
					end
				end
				function num=IndexStrArray_convert2Num(strArray,ind)
					num=NaN;
					try
						num=str2double(strArray(ind));
					end
				end
				function checkLength(a,b)
					if length(a)~=length(b)
						error("Length does not match")
					end
				end
			end
		end
	end
	%%
	properties
		DataStorageByCategory
	end
	methods
		function loadAllDataFromCategories(obj, options)
			arguments
				obj
				options.Properties2Store
				options.DisplayLoopInterval = []
				options.UseParfor (1,1) logical = false
			end
			SizeOfStructArray = size(obj.FilePathsList);
			ListOfFields = options.Properties2Store;
			data = InitializeStruct(SizeOfStructArray, ListOfFields);

			if options.UseParfor
				% Use parfor for parallel loading. Note: progress display is omitted.
				parfor i = 1:numel(obj.FilePathsList)
					if isempty(obj.FilePathsList(i))
						continue;
					end
					temp = load(obj.FilePathsList(i));
					for j = ListOfFields(:)'
						data(i).(j) = temp.obj.(j);
					end
				end
			else
				loopcounter = 0;
				for i = 1:numel(obj.FilePathsList)
					tic
					if isempty(obj.FilePathsList(i))
						continue;
					end
					temp = load(obj.FilePathsList(i));
					for j = ListOfFields(:)'
						data(i).(j) = temp.obj.(j);
					end
					loopcounter = loopcounter + 1;
					if isempty(options.DisplayLoopInterval)
						continue;
					end
					if loopcounter/options.DisplayLoopInterval == round(loopcounter/options.DisplayLoopInterval)
						fprintf("%d / %d ", loopcounter, numel(obj.FilePathsList));
						toc
					end
				end
			end
			obj.DataStorageByCategory = data;
			function data = InitializeStruct(SizeOfStructArray, ListOfFields)
				data = struct;
				ind = num2cell(SizeOfStructArray);
				for k = ListOfFields(:)'
					data(ind{:}).(k) = [];
				end
			end
		end
	end
	%% change trial start / end
	properties
		tstart_shift_sec
		tend_shift_sec
	end
	methods
		function modify_trialTime(obj,options)
			arguments
				obj
				options.tstart_shift_sec=0
				options.tend_shift_sec=0
			end
			if isempty(obj.cellac)|isempty(obj.DLC_raw)|isempty(obj.org)
				error("Execute CSP.loadData before this function")
			end
			obj.tstart_shift_sec=options.tstart_shift_sec;
			obj.tend_shift_sec=options.tend_shift_sec;
			ts=round(options.tstart_shift_sec*1000);
			te=round(options.tend_shift_sec*1000);
			neworg=obj.org;
			neworg.t_strial=neworg.t_strial+ts;
			neworg.t_etrial=neworg.t_etrial+te;
			validtrials=((neworg.t_strial>=1)&(neworg.t_etrial<=obj.org.t_etrial(end)));
			neworg.t_strial=neworg.t_strial(validtrials);
			neworg.t_etrial=neworg.t_etrial(validtrials);
			neworg.strial_t=[];
			neworg.etrial_t=[];
			neworg.motionLog=neworg.motionLog(validtrials);
			obj.loadData( ...
				UseCachedData					=true, ...
				Cached_cellac					=obj.cellac,...
				Cached_DLC_raw					=obj.DLC_raw,...
				Cached_org						=neworg, ...
				SamplingRate_Hz					=obj.SamplingRate_Hz, ...
				PeriodCutoff4Derivatives_sec	=obj.PeriodCutoff4Derivatives_sec);
		end
	end
	%% split tck to regions
	methods
		function all_ind=extract_cellINDs_of_regions(obj,RegNames)
			csP=copy(obj);
			Reg_log=csP.Reg_inRegLog_reg(containsSubstring(csP.Reg_RegionNames,RegNames));
			if sum(containsSubstring(csP.Reg_RegionNames,RegNames))~=numel(RegNames)
				warning("Region count mismatch")
			end
			Reg_log=cell2mat(Reg_log(:)');
			Reg_log=any(Reg_log,2);
			all_log=Reg_log(:)';
			all_ind=find(all_log);
		end
	end
	%% make time shifted data
	properties
		TimeShiftInSec
	end
	methods
		function makeTimeShiftedData(obj,options)
			arguments
				obj
				options.TimeShiftInSec	(1,1)=0
			end
			if isempty(obj.cellac)|isempty(obj.DLC_raw)|isempty(obj.org)
				error("Execute CSP.loadData before this function")
			end
			timeshiftms=round(1000*options.TimeShiftInSec);
			shifted_org=makeTimeShiftedOrg(obj.org,timeshiftms);
			obj.TimeShiftInSec=options.TimeShiftInSec;
			obj.loadData( ...
				UseCachedData					=true, ...
				Cached_org						=shifted_org, ...
				SamplingRate_Hz					=obj.SamplingRate_Hz, ...
				PeriodCutoff4Derivatives_sec	=obj.PeriodCutoff4Derivatives_sec);
			
			function neworg=makeTimeShiftedOrg(org,timeshiftms)
				[otherind,f2p_t_ind]=ShiftAndGetOverlappingVec([1:length(org.f2p_t)],[1:length(org.f2p_t)],timeshiftms);
				neworg=org;
				neworg.f2p_t=org.f2p_t(f2p_t_ind);	neworg.t_f2p=[];	neworg.trig2p=[];
				neworg.sp_t=org.sp_t(otherind,:);
				neworg.fdlc_t=org.fdlc_t(otherind);	neworg.t_fdlc=[];	neworg.trigfdlc=[];
				neworg=shift_t_trial_And_extract_valid_trials(neworg,otherind);
				function neworg=shift_t_trial_And_extract_valid_trials(neworg,t_ind)
					neworg.t_strial=neworg.t_strial-t_ind(1)+1;
					neworg.t_etrial=neworg.t_etrial-t_ind(1)+1;
					firstTrial=find(neworg.t_strial>=1,1,"first");
					lastTrial=find(neworg.t_etrial<=t_ind(end),1,"last");
					neworg.t_strial=neworg.t_strial([firstTrial:lastTrial]);
					neworg.t_etrial=neworg.t_etrial([firstTrial:lastTrial]);
					neworg.strial_t=[];
					neworg.etrial_t=[];
					neworg.motionLog=neworg.motionLog([firstTrial:lastTrial]);
				end
			end
		end
		%%
		function plotTimeShiftedCellacAndOthers(obj,options)
			arguments
				obj
				options.TimeShiftInSec_vec
				options.CellInd2plot (1,:) =[]
				options.SPparams2plot=[]
				options.motion2plot (1,1) string {mustBeMember(options.motion2plot, ...
					["rot", "tilt", "yaw", "zvert", ...
					"p8_dir1", "p8_dir2", "p8_dir3", ...
					"p8_dir4", "p8_dir5", "p8_dir6", ...
					"p8_dir7", "p8_dir8"])}="p8_dir8"
				options.BodyPart2plot (1,1) string {mustBeMember(options.BodyPart2plot, ...
					["RightHand", "RightElbow", "RightToe", ...
					"RightHeel", "RightKnee", "BackU", ...
					"BackM", "BackL", "Tail", "LeftHand", ...
					"LeftElbow", "LeftToe", "LeftHeel", ...
					"LeftKnee"])}="RightToe"
			end
			csp=copy(obj);
			timeshift=options.TimeShiftInSec_vec;
			for i=1:length(timeshift)
				csp.makeTimeShiftedData(TimeShiftInSec=timeshift(i));
				cellac=csp.cellac_mtype{csp.org.ind_mtypes(options.motion2plot)};
				SPparams=csp.sp_mtype{csp.org.ind_mtypes(options.motion2plot)};
				DLCparams=csp.DLC_mtype{csp.org.ind_mtypes(options.motion2plot)};
				DLCparams=DLCparams.tXYZbp{csp.rec3d.ind_bp(options.BodyPart2plot)};
				nexttile;
				cellind2plot=options.CellInd2plot;
				if isempty(cellind2plot)
					cellind2plot=1;
				end
				SPparams2plot=options.SPparams2plot;
				if isempty(SPparams2plot)
					SPparams2plot=[1:2];
				end
				DLCcoord2plot=1:3;
				spparamlist=["theta","phi","tau","x","y","z"];
				DLClist=["x","y","z"];
				% plotMultiDim([cellac(:,cellind2plot,:),SPparams(:,SPparams2plot,:),DLCparams(:,DLCcoord2plot,:)]);
				plotSimilarActivity([cellac(:,cellind2plot,:),SPparams(:,SPparams2plot,:),DLCparams(:,DLCcoord2plot,:)], labels=[arrayfun(@(a)sprintf("Cell activity of cell #%d",a),cellind2plot), ...
					arrayfun(@(a)sprintf("SPparams: %s",spparamlist(a)),SPparams2plot), ...
					arrayfun(@(a)sprintf("%s coordinate of %s",DLClist(a),options.BodyPart2plot),DLCcoord2plot)], ...
					normalizeData=true, ...
					offset_prctile=100);
				% legend([arrayfun(@(a)sprintf("Cell activity of cell #%d",a),cellind2plot), ...
				% 	arrayfun(@(a)sprintf("SPparams: %s",spparamlist(a)),SPparams2plot), ...
				% 	arrayfun(@(a)sprintf("%s coordinate of %s",DLClist(a),options.BodyPart2plot),DLCcoord2plot)]);
				title(sprintf("Time shift by %g sec",timeshift(i)))
				grid on;
				drawnow;
			end
			ChangeFigureSettings(MasterTitle=sprintf("Motion: %s",options.motion2plot),MasterTitleSize=20,TitleSize=15,Linkaxes="xy",LegendLoc="east",Legend=[arrayfun(@(a)sprintf("Cell activity of cell #%d",a),cellind2plot), ...
					arrayfun(@(a)sprintf("SPparams: %s",spparamlist(a)),SPparams2plot), ...
					arrayfun(@(a)sprintf("%s coordinate of %s",DLClist(a),options.BodyPart2plot),DLCcoord2plot)], ...
					LegendSize=13)
		end
		
		%% load_specified
		function load_specified(obj, plane, option)
			% LOAD Loads the specified property/properties for a given plane.
			%
			% Usage:
			%   obj.load(plane, "p8_tckd", "Skew")
			%
			% Valid property names (as strings) are:
			%   "p8_tckd", "rot_tck", "tilt_tck", "yaw_tck", "zvert_tck", ...
			%   "p8_tpkd", "rot_tpk", "tilt_tpk", "yaw_tpk", "zvert_tpk", ...
			%   "p8_DLC", "rot_DLC", "tilt_DLC", "yaw_DLC", "zvert_DLC", ...
			%   "rec3d", "Skew", "regROI"
			arguments
				obj
				plane
			end
			arguments (Repeating)
				option (1,:) string {mustBeMember(option, ["p8_tckd","rot_tck","tilt_tck","yaw_tck","zvert_tck", ...
					"p8_tpkd","rot_tpk","tilt_tpk","yaw_tpk","zvert_tpk", ...
					"p8_DLC","rot_DLC","tilt_DLC","yaw_DLC","zvert_DLC", ...
					"rec3d","Skew","regROI"])}
			end
			obj.plane=plane;
			% Loop over each requested property.
			for k = 1:length(option)
				prop = option{k};
				switch prop
					% case "p8_tckd"
					% 	temp = load(obj.info.NoAltplane.dir_p8(plane));
					% 	obj.p8_tckd = temp.p8_tckd;
					% case "rot_tck"
					% 	temp = load(obj.info.NoAltplane.dir_rot(plane));
					% 	obj.rot_tck = temp.rot_tck;
					% case "tilt_tck"
					% 	temp = load(obj.info.NoAltplane.dir_tilt(plane));
					% 	obj.tilt_tck = temp.tilt_tck;
					% case "yaw_tck"
					% 	temp = load(obj.info.NoAltplane.dir_yaw(plane));
					% 	obj.yaw_tck = temp.yaw_tck;
					% case "zvert_tck"
					% 	temp = load(obj.info.NoAltplane.dir_zvert(plane));
					% 	obj.zvert_tck = temp.zvert_tck;
					% case "p8_tpkd"
					% 	temp = load(obj.info.NoAltplane.dir_SPp8(plane));
					% 	obj.p8_tpkd = temp.p8_tpkd;
					% case "rot_tpk"
					% 	temp = load(obj.info.NoAltplane.dir_SProt(plane));
					% 	obj.rot_tpk = temp.rot_tpk;
					% case "tilt_tpk"
					% 	temp = load(obj.info.NoAltplane.dir_SPtilt(plane));
					% 	obj.tilt_tpk = temp.tilt_tpk;
					% case "yaw_tpk"
					% 	temp = load(obj.info.NoAltplane.dir_SPyaw(plane));
					% 	obj.yaw_tpk = temp.yaw_tpk;
					% case "zvert_tpk"
					% 	temp = load(obj.info.NoAltplane.dir_SPzvert(plane));
					% 	obj.zvert_tpk = temp.zvert_tpk;
					% case "p8_DLC"
					% 	temp = load(obj.info.NoAltplane.dir_DLCp8(plane));
					% 	obj.p8_DLC = temp.DLC_p8;
					% case "rot_DLC"
					% 	temp = load(obj.info.NoAltplane.dir_DLCrot(plane));
					% 	obj.rot_DLC = temp.DLC_rot;
					% case "tilt_DLC"
					% 	temp = load(obj.info.NoAltplane.dir_DLCtilt(plane));
					% 	obj.tilt_DLC = temp.DLC_tilt;
					% case "yaw_DLC"
					% 	temp = load(obj.info.NoAltplane.dir_DLCyaw(plane));
					% 	obj.yaw_DLC = temp.DLC_yaw;
					% case "zvert_DLC"
					% 	temp = load(obj.info.NoAltplane.dir_DLCzvert(plane));
					% 	obj.zvert_DLC = temp.DLC_zvert;

					case "rec3d"
						temp = load(fullfile(obj.dirlist.DLCdata_info(plane)));			rec3d = temp.rec3d;
						rec3d.ind_bp=@(b)arrayfun(@(a)find(rec3d.bp_ind==a),b);
						rec3d.name_ind=rec3d.bp_ind;
						rec3d.ind_name=rec3d.ind_bp;
						rec3d.Dotname_ind=rec3d.jointList;
						obj.rec3d = rec3d;
					case "Skew"
						temp = load(fullfile(obj.dirlist.skewness(plane)));				skew  = temp.skew;
						obj.Skew=skew;
					case "regROI"
						temp = load(fullfile(obj.dirlist.reg(plane)));					reg   = temp.reg;
						obj.regROI=reg.roi;
						obj.calculatemidline();
						obj.makeHalfAtlas();
						obj.makeHalfCoord();
					otherwise
						obj.loadData(plane);
				end
			end

			% Optionally, update cell activity normalization:
			% obj.normalizeCellac();

			% Preprocess the SP data.
			try
				for j=["p8_tpkd","rot_tpk","tilt_tpk","yaw_tpk","zvert_tpk"]
					obj.(j) = PreprocessSPdata(obj.(j));
				end
			catch
				% warning("SP data preprocessing skipped")
			end

			% assess hemisphere
			try
				obj.assessInsideBrainRegions();
			catch
				% warning("BrainRegion assessment skipped")
			end
		end


		function normalizeCellac(obj)
			% Normalize cell data
			allcellac = {obj.p8_tckd; obj.rot_tck; obj.tilt_tck; obj.yaw_tck; obj.zvert_tck};
			% Concatenate along dimensions 3 to ndims(a) and reduce to an n x c matrix
			allcellac = cellfun(@(a) multiBlockCat(a, 3:ndims(a), 1), allcellac, 'UniformOutput', false);
			allcellac = cell2mat(allcellac);
			cell_sigma = std(allcellac, 0, 1);
			cell_mu    = mean(allcellac, 1);
			f = @(a) (a - cell_mu) ./ cell_sigma;
            
            % Normalize each cell data property
            obj.p8_tckd   = f(obj.p8_tckd);
            obj.rot_tck   = f(obj.rot_tck);
            obj.tilt_tck  = f(obj.tilt_tck);
			obj.yaw_tck   = f(obj.yaw_tck);
			obj.zvert_tck = f(obj.zvert_tck);
		end
		%%
		function generatePlaneIND4eachMouse(obj)
			uniqueStr=unique(obj.info.NoAltplane.nameMouse_p,'stable');
			obj.plane_MouseIND=[];
			obj.mouseName_MouseIND=string;
			for i=1:length(uniqueStr)
				ind=find(obj.info.NoAltplane.nameMouse_p==uniqueStr(i));
				obj.plane_MouseIND{i}=ind;
				obj.mouseName_MouseIND(i)=uniqueStr(i);
			end
		end
		%% assess Brain regions %%%%%%%%%%%%%%%
		function assessInsideBrainRegions(obj)
			str=assessLeftOrRight();
			obj.Reg_LeftOrRight_reg=str;
			if str=="Left"
				obj.Reg_contour_reg=obj.contour_LeftRegion;
			elseif str=="Right"
				obj.Reg_contour_reg=obj.contour_RightRegion;
			end

			ROI=obj.regROI;
			Contour_g=obj.Reg_contour_reg;
			inRegion_g=assessInsideRegion_g(ROI,Contour_g);
			obj.Reg_inRegLog_reg=inRegion_g;
			obj.Reg_RegionNames=convertCharsToStrings(obj.RegionNames_type);
			function str=assessLeftOrRight()
				if obj.info.NoAltplane.logLeft_p(obj.plane)==1 && obj.info.NoAltplane.logRight_p(obj.plane)==0
					str="Left";
				elseif obj.info.NoAltplane.logLeft_p(obj.plane)==0 && obj.info.NoAltplane.logRight_p(obj.plane)==1
					str="Right";
				else
					error("Cannot assign hemispheres")
				end
			end
		end
		%%
		function plotBrainRegion(obj,options)
			arguments
				obj
				options.BrainRegion2plot (1,1) string {mustBeMember(options.BrainRegion2plot, ...
            ["PPC", "M2", "M1", "Trunk", "Whisker", "V1", "UL", "LL"])}
			end
			inRegion_g=obj.Reg_inRegLog_reg;
			a=convertCharsToStrings(obj.Reg_RegionNames);
			contourind=find(a==options.BrainRegion2plot);
			ROI=obj.regROI(inRegion_g{contourind},:);
			plotContour(roi=ROI,contour=obj.Reg_contour_reg{contourind},img=obj.AtlasImage);
		end

		%%
		function obj = clearData(obj,options)
			arguments
				obj
				options.Properties2Erase (1,:) string = []
				options.Properties2keep (1,:) string = []
			end

			% Get all non-constant properties
			mc = metaclass(obj);
			allProps = {mc.PropertyList.Name};
			allProps([mc.PropertyList.Constant]) = [];  % remove constants

			% Properties that must never be cleared
			defaultKeep = ["plane","dirlist","info","dir_BrainRegions_used","dir_AllMiceDataFile_used"];

			% Determine which properties to reset based on inputs:
			% - If Properties2keep is given, keep those (plus defaultKeep), clear the rest.
			% - Else if Properties2Erase is given, clear only those not in defaultKeep.
			% - Otherwise, clear all properties except defaultKeep.
			if ~isempty(options.Properties2keep)
				toErase = setdiff(allProps, [options.Properties2keep, defaultKeep]);
			elseif ~isempty(options.Properties2Erase)
				toErase = setdiff(options.Properties2Erase, defaultKeep);
			else
				toErase = setdiff(allProps, defaultKeep);
			end

			% Create a new object to retrieve the default property values
			% tempObj = feval(class(obj));
			tempObj=CSP(dir_AllMiceDataFile_used=obj.dir_AllMiceDataFile_used, ...
				dir_BrainRegions_used=obj.dir_BrainRegions_used);

			% Reset each property to its default value
			for p = toErase
				if isprop(obj, p)
					obj.(p) = tempObj.(p);
				end
			end
		end
		%%
		function clearAndSave(obj,options)
			arguments
				obj
				options.Properties2Erase (1,:) string = []
				options.Properties2keep (1,:) string = []
				options.filename =[]
				options.additionalVariable2save=[];
				options.SaveAs_v7p3 (1,1) logical =false
			end
			additionalVariable2save=options.additionalVariable2save;
			filename=options.filename;
			SaveAs_v7p3=options.SaveAs_v7p3;
			options=rmfield(options,["filename","additionalVariable2save","SaveAs_v7p3"]);
			options=namedargs2cell(options);
			obj.clearData(options{:});
			if isempty(additionalVariable2save)
				if SaveAs_v7p3
					save(filename,"obj","-v7.3");
				else
					save(filename,"obj");
				end
			else
				if SaveAs_v7p3
					save(filename,"obj","additionalVariable2save","-v7.3");
				else
					save(filename,"obj","additionalVariable2save");
				end
			end
		end



		%%
		function saveObj(obj, filename)
			if ~exist(foldername, 'dir')
				mkdir(foldername);
			end

			% Create a fresh object to get each property’s default
			% tempObj = feval(class(obj));
			tempObj=CSP(dir_AllMiceDataFile_used=obj.dir_AllMiceDataFile_used, ...
				dir_BrainRegions_used=obj.dir_BrainRegions_used);

			savedProps = {};
			unsavedProps = {};

			propNames = properties(obj);
			for i = 1:numel(propNames)
				propName = propNames{i};
				try
					value      = obj.(propName);
					defaultVal = tempObj.(propName);

					% Skip if exactly equal to default
					if isequal(value, defaultVal)
						continue
					end

					% Save non‑default value
					eval([propName ' = value;']);
					filename = fullfile(foldername, propName + ".mat");
					save(filename, propName);
					clear(propName)

					savedProps{end+1} = propName;
				catch ME
					warning('Failed to save "%s": %s', propName, ME.message);
					unsavedProps{end+1} = propName;
				end
			end

			metadata.saved   = savedProps;
			metadata.unsaved = unsavedProps;
			save(fullfile(foldername, "metadata.mat"), 'metadata');

			disp('Saved properties:'), disp(savedProps)
			disp('Failed to save:'), disp(unsavedProps)
		end
	end
	methods
		function ntl=check_and_assign_nested_tiledlayout(obj,ntl,options)
			arguments
				obj
				ntl
				options.make_new_if_empty (1,1) logical =true
			end
			if isempty(ntl)
				if options.make_new_if_empty
					warning("No nested_tiledlayout detected. creating a new one...")
					ntl=nested_tiledlayout();
				else
					error("Specify nested_tiledlayout object as an input");
				end
			end
		end
	end

%% Define the current data 2 use
	properties
		tck
		tpk
		DLC

		SP_T		% transformation matrix of current SP params
		invDLC	% DLC but XYZ is made so that the ground is flat
		DLC_allNumeric		% fields with cell --> cell is converted to mat
		invDLC_allNumeric	% fields with cell --> cell is converted to mat
	end

	methods
		%%
        function setCurrentData(obj, option)
            arguments
                obj
                option (1,1) string {mustBeMember(option, ["p8","rot","tilt","yaw","zvert"])}
            end
            switch option
                case "p8"
                    obj.tck = obj.p8_tckd;
                    obj.tpk = obj.p8_tpkd;
                    obj.DLC = obj.p8_DLC;
                case "rot"
                    obj.tck = obj.rot_tck;
                    obj.tpk = obj.rot_tpk;
                    obj.DLC = obj.rot_DLC;
                case "tilt"
                    obj.tck = obj.tilt_tck;
                    obj.tpk = obj.tilt_tpk;
                    obj.DLC = obj.tilt_DLC;
                case "yaw"
                    obj.tck = obj.yaw_tck;
                    obj.tpk = obj.yaw_tpk;
                    obj.DLC = obj.yaw_DLC;
                case "zvert"
                    obj.tck = obj.zvert_tck;
                    obj.tpk = obj.zvert_tpk;
                    obj.DLC = obj.zvert_DLC;
			end

			obj.make_DLC_speed();
			obj.make_DLC_accel();

			% Preprocess the SP data.
			obj.tpk = PreprocessSPdata(obj.tpk);

			% obtain T
			obj.SP_T=SPparams2TransMat(obj.tpk);

			% calculate DLC inverse transformed by stewart platform
			% obj.invDLC = applyInverseSPTrans2DLC(obj.DLC, obj.SP_T,1/obj.SamplingRate_Hz);

			% obtian all numeric DLC data (still a struct though)
			obj.DLC_allNumeric=obj.DLC;
			DLCflds=string(fields(obj.DLC_allNumeric));
			for k=DLCflds(:)'
				try
					obj.DLC_allNumeric.(k)=cell2mat(obj.DLC_allNumeric.(k));
				end
			end

			% 2026/04/22
			% % % % % obtian all numeric invDLC data (still a struct though)
			% % % % obj.invDLC_allNumeric=obj.invDLC;
			% % % % invDLCflds=string(fields(obj.invDLC_allNumeric));
			% % % % for k=invDLCflds(:)'
			% % % % 	try
			% % % % 		obj.invDLC_allNumeric.(k)=cell2mat(obj.invDLC_allNumeric.(k));
			% % % % 	end
			% % % % end	% concatenate data

			obj.DLC_inv_or_normal="normalPosture";
		end
		%%
		function smoothDLC(obj,options)
			arguments
				obj
				options
			end
		end
		%%
		function convertedDLC=DLC_convert(obj,DLCdata2convert,option)
			arguments
				obj
				DLCdata2convert
				option (1,1) string {mustBeMember(option,["Default","AllNumeric","All2ndDimToCell_exceptXYZcoordMaintained","All2ndDImToCell"])}
			end
			% just first convert all to numeric first
			flds=string(fields(DLCdata2convert));
			for k=flds(:)'
				if iscell(DLCdata2convert.(k))
					DLCdata2convert.(k)=cell2mat(DLCdata2convert.(k)(:)');
				end
			end

			% if All2ndDImToCell was selected
			specialcase=false;
			if option=="All2ndDImToCell"
				option="AllNumeric";
				specialcase=true;
			end


			XYZ=["tXYZbp","tDXYZbp","tD2XYZbp"];
			Ang=["tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"];

			convertedDLC=[];
			flds=string(fields(DLCdata2convert));
			for k=flds(:)'
				if ismember(k,XYZ)
					currdata=DLCdata2convert.(k);
					newdata=[];
					if iscell(currdata)
						if length(currdata)~=length(obj.rec3d.name_ind)...14
							error('XYZ data needs 14 elements if it is a cell')
						end
						if option=="AllNumeric"
							newdata=cell2mat(currdata(:)');
						else
							newdata=currdata;
						end
					else
						if size(currdata,2)~=3*length(obj.rec3d.name_ind)...14
							error('XYZ data needs size of 42 in the 2nd dimension')
						end
						if option=="AllNumeric"
							newdata=currdata;
						else
							IND=repmat({':'},1,ndims(currdata));
							newdata=[];
							for i=1:length(obj.rec3d.name_ind)
								IND{2}=[3*(i-1)+1:3*i];
								newdata{i}=currdata(IND{:});
							end
						end
					end
					convertedDLC.(k)=newdata;
				elseif ismember(k,Ang)
					currdata=DLCdata2convert.(k);
					newdata=[];
					if iscell(currdata)
						if length(currdata)~=size(obj.rec3d.Dotname_ind,1)...16
							error('Ang data needs 16 elements if it is a cell')
						end
						if option~="All2ndDimToCell_exceptXYZcoordMaintained"
							newdata=cell2mat(currdata(:)');
						else
							newdata=currdata;
						end
					else
						if size(currdata,2)~=size(obj.rec3d.Dotname_ind,1)...16
							error('Ang data needs size of 16 in the 2nd dimension')
						end
						if option~="All2ndDimToCell_exceptXYZcoordMaintained"
							newdata=currdata;
						else
							IND=repmat({':'},1,ndims(currdata));
							newdata=[];
							for i=1:size(obj.rec3d.Dotname_ind,1)
								IND{2}=i;
								newdata{i}=currdata(IND{:});
							end
						end
					end
					convertedDLC.(k)=newdata;
				elseif k=="triError"
					currdata=DLCdata2convert.(k);
					newdata=[];
					if iscell(currdata)
						if length(currdata)~=length(obj.rec3d.name_ind)...14
							error('TriError needs 14 elements if it is a cell')
						end
						if option~="All2ndDimToCell_exceptXYZcoordMaintained"
							newdata=cell2mat(currdata(:)');
						else
							newdata=currdata;
						end
					else
						if size(currdata,2)~=length(obj.rec3d.name_ind)...14
							error('TriError needs size of 14 in the 2nd dimension')
						end
						if option~="All2ndDimToCell_exceptXYZcoordMaintained"
							newdata=currdata;
						else
							IND=repmat({':'},1,ndims(currdata));
							newdata=[];
							for i=1:length(obj.rec3d.name_ind)...14
								IND{2}=i;
								newdata{i}=currdata(IND{:});
							end
						end
					end
					convertedDLC.(k)=newdata;
				else
					error(append("The field ",k," is weird. What is this?"))
				end
			end

			if specialcase
				for k=flds(:)'
					currdata=convertedDLC.(k);
					IND=repmat({':'},1,ndims(currdata));
					convertedDLC.(k)=cell(1,size(currdata,2));
					for i=1:size(currdata,2)
						IND{2}=i;
						convertedDLC.(k){i}=currdata(IND{:});
					end
				end
			end

			%%%% use code below to make sure this is right
			% A=CSP(info,plane=1);
			% A.loadData();
			% A.setCurrentData("p8");
			% a=[];
			% a.Default=A.DLC_convert(A.DLC,"Default");
			% a.All2ndDimToCell_exceptXYZcoordMaintained=A.DLC_convert(A.DLC,"All2ndDimToCell_exceptXYZcoordMaintained");
			% a.AllNumeric=A.DLC_convert(A.DLC,"AllNumeric");
			% a.All2ndDImToCell=A.DLC_convert(A.DLC,"All2ndDImToCell");
			% currstr="All2ndDImToCell"
			% b=A.DLC_convert(b,currstr)
			% isequal(a.Default,A.DLC)
			% isequal(a.AllNumeric,A.DLC_allNumeric)
			% isequal(b,a.(currstr))
		end
		%%
		function DLC=DLC_convert2numeric(obj,DLC)
			DLC=structfun(@(a)custom_cell2mat(a),DLC,UniformOutput=false);
			function newarray=custom_cell2mat(array)
				try
					newarray=cell2mat(array);
				catch
					newarray=array;
				end
			end
		end
		function DLC=DLC_convert2cell(obj,numeric_DLC,options)
			arguments
				obj
				numeric_DLC
				options.convert_all=false
			end
			cellflds=["tXYZbp","tDXYZbp","tD2XYZbp"];
			flds=string(fields(numeric_DLC));
			DLC=numeric_DLC;
			for i=flds(:)'
				if any(containsSubstring(cellflds,i))
					data=numeric_DLC.(i);
					IND=num2cell(size(data));
					IND{2}=repmat(3,1,size(data,2)/3);
					newdata=mat2cell(data,IND{:});
					DLC.(i)=newdata;
				elseif options.convert_all
					data=numeric_DLC.(i);
					IND=num2cell(size(data));
					IND{2}=ones(1,size(data,2));
					newdata=mat2cell(data,IND{:});
					DLC.(i)=newdata;
				end
			end
		end
		function indexed_array=index_numericDLC(obj,numericDLC,fld,index_2nd_dim)
			cellflds=["tXYZbp","tDXYZbp","tD2XYZbp"];
			if any(containsSubstring(cellflds,fld))
				IND=3*(index_2nd_dim-1)+1:3*index_2nd_dim;
			else
				IND=index_2nd_dim;
			end
			allind=createIND(numericDLC.(fld),2,IND);
			indexed_array=numericDLC.(fld)(allind{:});
		end
	end
	methods (Static)
		function catdata=cat3Dcoordinates(data)
			if ~iscell(data)
				error("data has to be cell")
			end
			catdata=cat(2,data{:});
		end
		function splitdata=split3Dcoordinates(data)
			if ~isnumeric(data)
				error("data has to be a numeric array with 2nd dimension being concatenation of 3d coordinates")
			end
			splitCell=num2cell(size(data));
			if splitCell{2}/3-round(splitCell{2}/3)~=0
				error("2nd dimension has to be multiple of 3")
			end
			splitCell{2}=3+zeros(1,splitCell{2}/3);
			splitdata=mat2cell(data,splitCell{:});
		end
		%% Unit test
		% a=csp.cat3Dcoordinates(csp.DLC.tXYZbp);
		% b=csp.split3Dcoordinates(a);
		% isequal(b,csp.DLC.tXYZbp)
	end
	%% posture plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		rec3d_full
		bp2d
		camera_bp2d
		bp_bp2d
		bp_cam
		coord_cam
		vid_posture
	end
	methods
		function load_posture_video(obj,options)
			arguments
				obj
				options.rec3d=[]
			end
			filelist=listAllFiles(obj.dir_AllMiceDirectories);
			fileind=containsSubstring(filelist,obj.currPlaneName4Directory);
			temp=load(filelist(fileind));
			if isempty(options.rec3d)
				temp=load(temp.allDirectories.DLCdata);
				obj.rec3d_full=temp.rec3d;
			else
				obj.rec3d_full=options.rec3d;
			end
			rec3d=obj.rec3d_full;

			vid=VideoReader(rec3d.DLCvid);
			I=read(vid,1);
			bp2d = [
				"RightHand1",   "RightElbow1",  "RightToe1",    "RightHeel1",   "RightKnee1", ...
				"BackU1",       "BackM1",       "BackL1",       "Tail1",        "LeftHand2", ...
				"LeftElbow2",   "LeftToe2",     "LeftHeel2",    "LeftKnee2",    "BackU2", ...
				"BackM2",       "BackL2",       "Tail2",        "LeftHand3",    "LeftElbow3", ...
				"LeftToe3",     "LeftHeel3",    "LeftKnee3",    "RightHand4",   "RightElbow4", ...
				"RightToe4",    "RightHeel4",   "RightKnee4"
				];
			lastChars = extractAfter(bp2d, strlength(bp2d)-1);
			camera_bp2d = str2double( lastChars );
			bp_bp2d = extractBefore(bp2d, strlength(bp2d));

			bp_cam={};
			coord_cam={};
			for cam=1:4
				bp_cam{cam}=bp_bp2d(camera_bp2d==cam);
				coord_cam{cam}=cell(1,numel(bp_cam{cam}));
			end
			% get image coordinates for 2d bp
			for bp=1:numel(obj.rec3d.bp_ind)
				for i=1:numel(rec3d.imageCoord_bp{bp})
					cam=obj.assignCamera(rec3d.imageCoord_bp{bp}{i},size(I));
					ind=(bp_cam{cam}==obj.rec3d.bp_ind(bp));
					coord_cam{cam}{ind}=rec3d.imageCoord_bp{bp}{i};
				end
			end
			obj.bp2d=bp2d;
			obj.camera_bp2d=camera_bp2d;
			obj.bp_bp2d=bp_bp2d;
			obj.bp_cam=bp_cam;
			obj.coord_cam=coord_cam;
			obj.vid_posture=VideoReader(rec3d.DLCvid);
		end 
		function imshow_posture(obj,options)
			arguments
				obj
				options.IND
				options.imscaler=1
			end
			rec3d		  =obj.rec3d_full;
			bp2d          = obj.bp2d;
			camera_bp2d   = obj.camera_bp2d;
			bp_bp2d       = obj.bp_bp2d;
			bp_cam        = obj.bp_cam;
			coord_cam     = obj.coord_cam;
			timepoint=options.IND;
			immax=[200,200,200,97];
			immax=immax/255;
			orgtime=obj.timestamp(timepoint);
			frameIND=round(obj.org.fdlc_t(orgtime));
			vid=obj.vid_posture;
			I=read(vid,frameIND);
			I=imresize(I,options.imscaler);
			imshow(obj.adjustBrightness(I,immax));
		end
		function scatter_posture(obj,options)
			arguments
				obj
				options.IND
				options.dotsize=100
				options.plotNose=false
				options.imscaler=1
			end
			rec3d		  =obj.rec3d_full;
			bp2d          = obj.bp2d;
			camera_bp2d   = obj.camera_bp2d;
			bp_bp2d       = obj.bp_bp2d;
			bp_cam        = obj.bp_cam;
			coord_cam     = obj.coord_cam;
			timepoint=options.IND;
			dotsize=options.dotsize;
			orgtime=obj.timestamp(timepoint);
			frameIND=round(obj.org.fdlc_t(orgtime));
			if options.plotNose
				cmap=turbo(length(obj.rec3d.bp_ind)+1);
			else
				cmap=turbo(length(obj.rec3d.bp_ind));
			end
			% scatter
			x=[];
			y=[];
			color2d={};
			for bp2dIND=1:numel(bp2d)
				cam=camera_bp2d(bp2dIND);
				bp=bp_bp2d(bp2dIND);
				coord=coord_cam{cam}{bp_cam{cam}==bp}(frameIND,:)*options.imscaler;
				x(end+1)=coord(1);
				y(end+1)=coord(2);
				color2d{end+1}=cmap(obj.rec3d.ind_bp(bp),:);
				find(bp_cam{cam}==bp);
			end
			scatter(x,y,dotsize,cat(1,color2d{:}),"filled");
		end
		function plot_posture(obj,options)
			arguments
				obj
				options.IND
				options.edgewidth=3
				options.plotNose=false
				options.imscaler=1;
			end
			rec3d		  =obj.rec3d_full;
			bp2d          = obj.bp2d;
			camera_bp2d   = obj.camera_bp2d;
			bp_bp2d       = obj.bp_bp2d;
			bp_cam        = obj.bp_cam;
			coord_cam     = obj.coord_cam;
			timepoint=options.IND;
			edgewidth=options.edgewidth;
			orgtime=obj.timestamp(timepoint);
			frameIND=round(obj.org.fdlc_t(orgtime));

			if options.plotNose
				cmap=turbo(length(obj.rec3d.bp_ind)+1);
			else
				cmap=turbo(length(obj.rec3d.bp_ind));
			end

			holdstate=ishold(gca());
			% plot
			alledges=[obj.rec3d.edgeList;"RightElbow","RightKnee";"LeftElbow","LeftKnee"];
			for cam=1:4
				currbps=bp_cam{cam};
				for edgeIND=1:size(alledges,1)
					edge=alledges(edgeIND,:);
					ind=[find(currbps==edge(1)),find(currbps==edge(2))];
					if length(ind)~=2
						continue;
					end
					coord=coord_cam{cam}(ind);
					color=mean(cmap(obj.rec3d.ind_bp(edge),:),1);
					x=[coord{1}(frameIND,1),coord{2}(frameIND,1)]*options.imscaler;
					y=[coord{1}(frameIND,2),coord{2}(frameIND,2)]*options.imscaler;
					plot(x,y,color=color,LineWidth=edgewidth);
					hold on
				end
			end
			if ~holdstate
				hold off
			end
		end
		function camIdx = assignCamera(obj,points, imgSize)
			% ASSIGNCAMERA  Determine which of 4 quadrants a set of points belongs to
			%
			%   camIdx = ASSIGNCAMERA(points, imgSize) takes:
			%     - points: an N×2 array of [x, y] image coordinates
			%     - imgSize: a two-element vector [height, width]
			%
			%   It computes the median point of all rows in points, then
			%   splits the image plane in half horizontally and vertically,
			%   assigning:
			%     1 = top-right quadrant
			%     2 = top-left quadrant
			%     3 = bottom-left quadrant
			%     4 = bottom-right quadrant
			%
			%   Example:
			%     pts = [10 20; 300 50; 400 10; 50 200];
			%     imgSize = [480, 640];
			%     cam = assignCamera(pts, imgSize)   % returns e.g. 1,2,3 or 4

			% Compute median x and y:
			m = median(points, 1);   % [medianX, medianY]
			cx = imgSize(2) / 2;     % center x
			cy = imgSize(1) / 2;     % center y

			% Determine quadrant:
			if m(1) > cx         % right half
				if m(2) <= cy    % top half
					camIdx = 1;  % top-right
				else
					camIdx = 4;  % bottom-right
				end
			else                % left half
				if m(2) <= cy    % top half
					camIdx = 2;  % top-left
				else
					camIdx = 3;  % bottom-left
				end
			end
		end
		function modI=adjustBrightness(obj,I,immax)
			I=rgb2gray(I);
			modI=I;
			halfy=size(modI,1)/2;
			halfx=size(modI,2)/2;
			i=1;
			modI(halfy+1:end,1:halfx,:)=imadjust(I(halfy+1:end,1:halfx,:),[0,immax(i)]);

			i=2;
			modI(1:halfy,1:halfx,:)=imadjust(I(1:halfy,1:halfx,:),[0,immax(i)]);

			i=3;
			modI(1:halfy,halfx+1:end,:)=imadjust(I(1:halfy,halfx+1:end,:),[0,immax(i)]);

			i=4;
			modI(halfy+1:end,halfx+1:end,:)=imadjust(I(halfy+1:end,halfx+1:end,:),[0,immax(i)]);
		end
	end
	properties
		handles_scatter3_posture
		handles_plot3_posture
	end
	methods
		function hdls=scatter3_posture(obj,posture_data,options)
			arguments
				obj
				posture_data	% either numeric or cell
				options.IND		% 1st dimension index to plot the posture
				options.DotSize=100
				options.plotNose=false
				options.xlim=[-90 20]
				options.ylim=[-40 40]
				options.zlim=[-15 40]
				options.setLim (1,1) logical =true
				options.ax=[]	% current ax
				options.plot_handles=[]
			end
			if isempty(options.ax)
				options.ax=gca;
			end
			% axes(options.ax);
			hdls=gobjects(1,0);
			hdl_ind=1;
			% nosedata=[32.8339,5.8359,16.8415];
			nosedata=[18.6827,-0.6206,16.2661];
			holdlog=ishold;
			if isnumeric(posture_data)
				posture_data=obj.split3Dcoordinates(posture_data);
			end
			if options.plotNose
				bp_color=turbo(length(obj.rec3d.bp_ind)+1);
			else
				bp_color=turbo(length(obj.rec3d.bp_ind));
			end

			xyz=[];
			for bp=1:length(obj.rec3d.bp_ind)
				xyz=cat(1,xyz,posture_data{bp}(options.IND,:));
			end
			[hdls,hdl_ind]=scatter3d(xyz(:,1),xyz(:,2),xyz(:,3),options,bp_color(1:end-1,:),hdls,hdl_ind);
			% scatter3(options.ax,xyz(:,1),xyz(:,2),xyz(:,3),options.DotSize,bp_color(1:end-1,:),"filled")
			if options.plotNose
				hold on
				[hdls,hdl_ind]=scatter3d(nosedata(1),nosedata(2),nosedata(3),options,bp_color(end,:),hdls,hdl_ind);
				% scatter3(nosedata(1),nosedata(2),nosedata(3),options.DotSize,bp_color(end,:),"filled")
			end
			if options.setLim
				xlim(options.xlim);
				ylim(options.ylim);
				zlim(options.zlim);
			end
			daspect([1,1,1])
			modes = {'off','on'};     % modes{1} = 'off', modes{2} = 'on'
			hold( modes{ holdlog+1 } );  % ishold returns 0 or 1 → index 1 or 2
			if ~isempty(options.plot_handles)
				hdls=options.plot_handles;
			end

			function [hdls,hdl_ind]=scatter3d(x,y,z,options,col,hdls,hdl_ind)
				if isempty(options.plot_handles)
					hdls(hdl_ind)=scatter3(options.ax,x,y,z,options.DotSize,col,"filled");
				else
					set(options.plot_handles(hdl_ind),"XData",x);
					set(options.plot_handles(hdl_ind),"YData",y);
					set(options.plot_handles(hdl_ind),"ZData",z);
				end
				hdl_ind=hdl_ind+1;
			end
		end
		function hdls=scatter3_posture_allPoints(obj,posture_data,options)
			arguments
				obj
				posture_data	% either numeric or cell
				options.alpha=0.3
				options.DotSize=1
				options.plotNose=false
				options.xlim=[-90 20]
				options.ylim=[-40 40]
				options.zlim=[-15 40]
				options.setLim (1,1) logical =true
				options.ax=[]	% current ax
				options.plot_handles=[]
			end
			if isempty(options.ax)
				options.ax=gca;
			end
			% axes(options.ax);
			hdls=gobjects(1,0);
			hdl_ind=1;
			% nosedata=[32.8339,5.8359,16.8415];
			nosedata=[18.6827,-0.6206,16.2661];
			holdlog=ishold;
			if isnumeric(posture_data)
				posture_data=obj.split3Dcoordinates(posture_data);
			end
			if options.plotNose
				bp_color=turbo(length(obj.rec3d.bp_ind)+1);
			else
				bp_color=turbo(length(obj.rec3d.bp_ind));
			end

			xyz=[];
			X=[];
			Y=[];
			Z=[];
			for bp=1:length(obj.rec3d.bp_ind)
				X=cat(2,X,posture_data{bp}(:,1));
				Y=cat(2,Y,posture_data{bp}(:,2));
				Z=cat(2,Z,posture_data{bp}(:,3));
				% xyz=cat(1,xyz,posture_data{bp}(options.IND,:));
			end
			[hdls,hdl_ind]=scatter3d(X,Y,Z,options,bp_color(1:end-1,:),hdls,hdl_ind);
			% scatter3(options.ax,xyz(:,1),xyz(:,2),xyz(:,3),options.DotSize,bp_color(1:end-1,:),"filled")
			if options.plotNose
				hold on
				[hdls,hdl_ind]=scatter3d(nosedata(1),nosedata(2),nosedata(3),options,bp_color(end,:),hdls,hdl_ind);
				% scatter3(nosedata(1),nosedata(2),nosedata(3),options.DotSize,bp_color(end,:),"filled")
			end
			if options.setLim
				xlim(options.xlim);
				ylim(options.ylim);
				zlim(options.zlim);
			end
			daspect([1,1,1])
			modes = {'off','on'};     % modes{1} = 'off', modes{2} = 'on'
			hold( modes{ holdlog+1 } );  % ishold returns 0 or 1 → index 1 or 2
			if ~isempty(options.plot_handles)
				hdls=options.plot_handles;
			end

			function [hdls,hdl_ind]=scatter3d(x,y,z,options,col,hdls,hdl_ind)
				if isempty(options.plot_handles)
					h=scatter3(options.ax,x,y,z,options.DotSize,col,"filled", ...
						"MarkerFaceAlpha",options.alpha,"MarkerEdgeAlpha",options.alpha);
					for i=1:length(h)
						hdls(hdl_ind)=h(i);
						hdl_ind=hdl_ind+1;
					end
				else
					set(options.plot_handles(hdl_ind),"XData",x);
					set(options.plot_handles(hdl_ind),"YData",y);
					set(options.plot_handles(hdl_ind),"ZData",z);
				end
				hdl_ind=hdl_ind+1;
			end
		end
		function hdls=plot3_posture(obj,posture_data,options)
			arguments
				obj
				posture_data	% either numeric or cell
				options.IND		% 1st dimension index to plot the posture
				options.edgesize=3
				options.plotNose=false
				options.xlim=[-90 20]
				options.ylim=[-40 40]
				options.zlim=[-15 40]
				options.setLim (1,1) logical =true
				options.ax=[]	% current ax
				options.plot_handles=[]
			end
			if isempty(options.ax)
				options.ax=gca;
			end
			% axes(options.ax);
			hdls=gobjects(1,0);
			hdl_ind=1;
			% nosedata=[32.8339,5.8359,16.8415];
			nosedata=[18.6827,-0.6206,16.2661];
			holdlog=ishold;
			if isnumeric(posture_data)
				posture_data=obj.split3Dcoordinates(posture_data);
			end
			if options.plotNose
				bp_color=turbo(length(obj.rec3d.bp_ind)+1);
			else
				bp_color=turbo(length(obj.rec3d.bp_ind));
			end
			edgecolor=[];
			for i=1:size(obj.rec3d.edgeList,1)
				bp1=obj.rec3d.ind_name(obj.rec3d.edgeList(i,1));
				bp2=obj.rec3d.ind_name(obj.rec3d.edgeList(i,2));
				edgecolor(i,:)=1/2*(bp_color(bp1,:)+bp_color(bp2,:));
			end
			additionalEdges=["RightElbow","RightKnee";"LeftElbow","LeftKnee"];
			additionaledgecolor=[];
			for i=1:size(additionalEdges,1)
				bp1=obj.rec3d.ind_name(additionalEdges(i,1));
				bp2=obj.rec3d.ind_name(additionalEdges(i,2));
				additionaledgecolor(i,:)=1/2*(bp_color(bp1,:)+bp_color(bp2,:));
			end
			nose_edges=[
				"BackU",      "nose";
				"RightElbow", "nose";
				"LeftElbow",  "nose"
				];
			nose_edgecolor=[];
			for i=1:size(nose_edges,1)
				bp1=obj.rec3d.ind_name(nose_edges(i,1));
				nose_edgecolor(i,:)=1/2*(bp_color(bp1,:)+bp_color(end,:));
			end

			X=[];
			Y=[];
			Z=[];

			for edge=1:size(obj.rec3d.edgeList,1)
				inds=obj.rec3d.ind_bp(obj.rec3d.edgeList(edge,:));
				bp1=posture_data{inds(1)}(options.IND,:);
				bp2=posture_data{inds(2)}(options.IND,:);
				[hdls,hdl_ind]=plot3d(bp1,bp2,options,edgecolor(edge,:),hdls,hdl_ind);
				% plot3(options.ax,[bp1(1);bp2(1)],[bp1(2);bp2(2)],[bp1(3);bp2(3)],"Color",edgecolor(edge,:),"LineWidth",options.edgesize);
				hold on;
			end
			for edge=1:size(additionalEdges,1)
				inds=obj.rec3d.ind_bp(additionalEdges(edge,:));
				bp1=posture_data{inds(1)}(options.IND,:);
				bp2=posture_data{inds(2)}(options.IND,:);
				[hdls,hdl_ind]=plot3d(bp1,bp2,options,additionaledgecolor(edge,:),hdls,hdl_ind);
				% plot3(options.ax,[bp1(1);bp2(1)],[bp1(2);bp2(2)],[bp1(3);bp2(3)],"Color",additionaledgecolor(edge,:),"LineWidth",options.edgesize);
				hold on;
			end
			if options.plotNose
				for edge=1:size(nose_edges,1)
					inds=obj.rec3d.ind_bp(nose_edges(edge,1));
					bp1=posture_data{inds(1)}(options.IND,:);
					bp2=nosedata;
					[hdls,hdl_ind]=plot3d(bp1,bp2,options,nose_edgecolor(edge,:),hdls,hdl_ind);
					% plot3(options.ax,[bp1(1);bp2(1)],[bp1(2);bp2(2)],[bp1(3);bp2(3)],"Color",nose_edgecolor(edge,:),"LineWidth",options.edgesize);
					hold on;
				end
			end
			modes = {'off','on'};     % modes{1} = 'off', modes{2} = 'on'
			hold( modes{ holdlog+1 } );  % ishold returns 0 or 1 → index 1 or 2
			if options.setLim
				xlim(options.xlim);
				ylim(options.ylim);
				zlim(options.zlim);
			end
			daspect([1,1,1])
			if ~isempty(options.plot_handles)
				hdls=options.plot_handles;
			end

			function [hdls,hdl_ind]=plot3d(bp1,bp2,options,col,hdls,hdl_ind)
				if isempty(options.plot_handles)
					hdls(hdl_ind)=plot3(options.ax,[bp1(1);bp2(1)],[bp1(2);bp2(2)],[bp1(3);bp2(3)],"Color",col,"LineWidth",options.edgesize);
				else
					set(options.plot_handles(hdl_ind),"XData",[bp1(1);bp2(1)]);
					set(options.plot_handles(hdl_ind),"YData",[bp1(2);bp2(2)]);
					set(options.plot_handles(hdl_ind),"ZData",[bp1(3);bp2(3)]);
				end
				hdl_ind=hdl_ind+1;
			end
		end
		function hdls=plot3_posture_pcaAxes(obj,options)
			arguments
				obj
				options.posture_mean (1,:)
				options.posture_axes (1,:)	% axis of pca of posture multiplied by its variance
				options.edgesize=1
				options.plotNose=false
				options.xlim=[-90 20]
				options.ylim=[-40 40]
				options.zlim=[-15 40]
				options.setLim (1,1) logical =true
				options.ax=[]	% current ax
				options.plot_handles=[]
				options.axes_color=""
			end
			if isempty(options.ax)
				options.ax=gca;
			end
			% axes(options.ax);
			hdls=gobjects(1,0);
			hdl_ind=1;
			holdlog=ishold;
			p1=options.posture_mean+options.posture_axes;
			p2=options.posture_mean-options.posture_axes;
			p1=obj.split3Dcoordinates(p1);
			p2=obj.split3Dcoordinates(p2);
			if options.plotNose
				bp_color=turbo(length(obj.rec3d.bp_ind)+1);
			else
				bp_color=turbo(length(obj.rec3d.bp_ind));
			end
			% i=1;
			% X=[cellfun(@(a)a(i),p1);cellfun(@(a)a(i),p2)];
			% i=2;
			% Y=[cellfun(@(a)a(i),p1);cellfun(@(a)a(i),p2)];
			% i=3;
			% Z=[cellfun(@(a)a(i),p1);cellfun(@(a)a(i),p2)];

			for bp=1:length(obj.rec3d.bp_ind)
				i=1;
				x=[p1{bp}(i);p2{bp}(i)];
				i=2;
				y=[p1{bp}(i);p2{bp}(i)];
				i=3;
				z=[p1{bp}(i);p2{bp}(i)];
				if options.axes_color==""
					[hdls,hdl_ind]=plot3d(x,y,z,options,options.axes_color,hdls,hdl_ind);
				else
					[hdls,hdl_ind]=plot3d(x,y,z,options,bp_color(bp,:),hdls,hdl_ind);
				end
				
				% plot3(options.ax,x,y,z,"Color",bp_color(bp,:),"LineWidth",options.edgesize)
			end
			
			modes = {'off','on'};     % modes{1} = 'off', modes{2} = 'on'
			hold( modes{ holdlog+1 } );  % ishold returns 0 or 1 → index 1 or 2
			if options.setLim
				xlim(options.xlim);
				ylim(options.ylim);
				zlim(options.zlim);
			end
			daspect([1,1,1])
			if ~isempty(options.plot_handles)
				hdls=options.plot_handles;
			end

			function [hdls,hdl_ind]=plot3d(x,y,z,options,col,hdls,hdl_ind)
				if isempty(options.plot_handles)
					hdls(hdl_ind)=plot3(options.ax,x,y,z,"Color",col,"LineWidth",options.edgesize);
				else
					set(options.plot_handles(hdl_ind),"XData",x);
					set(options.plot_handles(hdl_ind),"YData",y);
					set(options.plot_handles(hdl_ind),"ZData",z);
				end
				hdl_ind=hdl_ind+1;
			end
		end
		
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% cache all planes
	properties
		csp_plane={}
	end
	methods
		function cache_csp(obj,options)
			arguments
				obj
				options.plane=[]
				options.properties_to_store (1,:) string =[]
				options.properties_to_delete (1,:) string =[]
			end
			csp=copy(obj);
			csp.clearData(Properties2keep=options.properties_to_store,Properties2Erase=options.properties_to_delete);
			if isempty(obj.plane) && isempty(options.plane);error("Specify obj.plane or plane");end
			if isempty(options.plane)
				if isempty(obj.plane)
					error("Specify obj.plane or plane");
				else
					obj.csp_plane{obj.plane}=csp;
				end
			else
				obj.csp_plane{options.plane}=csp;
			end
		end
		function cache_all_planes(obj,options)
			arguments
				obj
				options.func_sequence
				options.properties_to_store (1,:) string =[]
				options.properties_to_delete (1,:) string =[]
				options.display_progress (1,1) logical =true
			end
			for plane=1:length(obj.PlaneName)
				csp=copy(obj);
				csp.plane=plane;
				apply_func_sequence(csp,options.func_sequence);
				if ~isempty(options.properties_to_delete)||~isempty(options.properties_to_store)
					csp.clearData(Properties2keep=options.properties_to_store,Properties2Erase=options.properties_to_delete);
				end
				obj.csp_plane{plane}=csp;
				if options.display_progress
					fprintf("Plane %d out of %d complete\n",plane,length(obj.PlaneName))
				end
			end
			function apply_func_sequence(csp,func)
				for i=1:length(func)
					f=func{i};
					f(csp);
				end
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% direction selectivity
	methods
		%% plot to see if cell activity to preferred direction is consistent
		function plotRobustnessOfSelectivityAcrossTrial(obj,options,plotoptions)
			arguments
				obj
				options.pre_margin_in_sec=0
				options.post_margin_in_sec=0
				options.divide_by_sum (1,1) logical =false
				options.use_margin_4_baseline (1,1) logical =false
				options.subtract_prctile =[]

				plotoptions.pre_margin_in_sec_4_display=0
				plotoptions.post_margin_in_sec_4_display=0
				plotoptions.skew_thresh=[]
				plotoptions.ntl=[]
				plotoptions.num_cell_to_plot=[]
				plotoptions.plot_colorbar (1,1) logical =false
			end
			options=namedargs2cell(options);
			obj.assess_direction_selectivity(options{:});
			skewlog=obj.Skew>=plotoptions.skew_thresh;
			cellind=find(skewlog);
			tckd=generate_tckd_4_plot(plotoptions.pre_margin_in_sec_4_display, ...
				plotoptions.post_margin_in_sec_4_display);
			tckd=tckd(:,skewlog,:,:);
			score=obj.direction_selectivity.score(skewlog);
			maxdir=obj.direction_selectivity.direction(skewlog);
			pval=obj.direction_selectivity.pval(skewlog);

			[~,sorted_ind]=sort(score,"descend");
			ntl=obj.check_and_assign_nested_tiledlayout(plotoptions.ntl,make_new_if_empty=true);
			for i=sorted_ind(1:plotoptions.num_cell_to_plot)'
				ntl.nexttile();
				currscore=score(i);
				currmaxdir=maxdir(i);
				currtckd=tckd(:,i,:,currmaxdir);
				kt=permute(currtckd,[3,1,2,4]);
				imagesc(kt,[-Inf,2]);
				title(sprintf("cell #%d, direction %d",cellind(i),maxdir(i)));
				if plotoptions.plot_colorbar
					colorbar;
				end
			end
			function tckd=generate_tckd_4_plot(pre_mrgn,post_mrgn)
				csp=copy(obj);
				csp.loadData_with_trial_margin(pre_margin_in_sec=pre_mrgn,post_margin_in_sec=post_mrgn);
				csp.setCurrentData("p8");
				tckd=csp.tck;
			end
		end
		%%
		function plotRobustnessOfSelectivityAcrossTrial_allPlanes(obj,options)
			arguments
				obj
				options.NumCells2plot
				options.TileHandle=[]
				options.SkewThresh=[];
				options.PvalThresh=[];
				options.trial_prctile_2_use=[0,100]
				options.sc_prctile=[10,90]
				options.cut_bottom_prctile (1,1) logical =true
				options.use_margin_4_baseline (1,1) logical = false
			end


			ind=[1:length(obj.info.NoAltplane.nameFull_p)];
			% ind=[1:2];
			obj.plotTile_allPlane(@plotfunc,ind,Exec_nexttile=false)
			ChangeFigureSettings(MasterTitle=sprintf("Robustness (x=time, y=trial), Cells with Skewness>=%g",options.SkewThresh), ...
				MasterTitleSize=20);
			function plotfunc(i)
				% csp=CSP(dir_AllMiceDataFile_used=obj.dir_AllMiceDataFile_used, ...
				% 	dir_BrainRegions_used=obj.dir_BrainRegions_used);
				% csp.loadData(i);
				% csp.setCurrentData("p8");
				csp=copy(obj.csp_plane{i});
				fig=gcf;
				t=tiledlayout(fig.Children,'flow','TileSpacing','tight','Padding','tight');
				t.Layout.Tile=i;
				t.Title.String=obj.info.NoAltplane.nameFull_p(i);
				opts=options;
				opts.TileHandle=t;
				opts=namedargs2cell(opts);
				csp.plotRobustnessOfSelectivityAcrossTrial(opts{:});
				drawnow;
				i
			end
		end
		%%
	end
	properties
		direction_selectivity	% struct
	end
	methods
		function assess_direction_selectivity(obj,options)
			arguments
				obj
				options.pre_margin_in_sec=0
				options.post_margin_in_sec=0
				options.divide_by_sum (1,1) logical =false
				options.use_margin_4_baseline (1,1) logical =false
				options.subtract_prctile =[]
			end
			if options.use_margin_4_baseline && ~isempty(options.subtract_prctile)
				warning("use_margin_4_baseline and subtract_prctile is incompatible")
			end
			obj.direction_selectivity=struct;
			tckd=obj.p8_tckd;
			pval=assess_ANOVA_pval(tckd);
			[score,maxdir,tckd,CD]=assess_direction_selectivity(obj);
			obj.direction_selectivity.pval=pval;
			obj.direction_selectivity.score=score;
			obj.direction_selectivity.direction=maxdir;
			obj.direction_selectivity.plane=obj.plane;
			obj.direction_selectivity.CD=CD;
			% obj.direction_selectivity.tckd=tckd;
			
			
			function [score,maxdir,tckd,CD]=assess_direction_selectivity(obj)
				tckd=obj.p8_tckd;
				CD=squeeze(mean(tckd,[1,3]));
				csp=copy(obj);
				pre_mrgn=options.pre_margin_in_sec;
				post_mrgn=options.post_margin_in_sec;
				if pre_mrgn~=0 || post_mrgn~=0
					csp.loadData_with_trial_margin(pre_margin_in_sec=pre_mrgn,post_margin_in_sec=post_mrgn);
					csp.setCurrentData("p8");
					mtckd=csp.tck;
					pre_mrgn=pre_mrgn*obj.SamplingRate_Hz;
					post_mrgn=post_mrgn*obj.SamplingRate_Hz;
					base_tckd=mtckd([1:pre_mrgn,size(mtckd,1)-post_mrgn+1:size(mtckd,1)],:,:,:);
					base_tckd=mean(base_tckd,1);
					tckd=tckd-base_tckd;
					ckd=permute(mean(tckd,1),[2,3,4,1]);
					ckd(ckd<0)=0;
					CD=permute(mean(ckd,2),[1,3,2]);
				end
				if options.use_margin_4_baseline
					csp=copy(obj);
					mrgnsec=1;
					csp.loadData_with_trial_margin(pre_margin_in_sec=mrgnsec, ...
						post_margin_in_sec=mrgnsec);
					csp.setCurrentData("p8");
					mtckd=csp.tck;
					mrgn_duration=mrgnsec*obj.SamplingRate_Hz;
					base_tckd=mtckd([1:mrgn_duration,size(mtckd,1)-mrgn_duration:size(mtckd,1)], ...
						:,:,:);
					base_tckd=mean(base_tckd,1);
					tckd=tckd-base_tckd;
					tckd(tckd<0)=0;
					CD=squeeze(mean(tckd,[1,3]));
				end
				if ~isempty(options.subtract_prctile)
					TKDc=reshape(permute(tckd,[2,1,3,4]),size(tckd,2),[])';
					plist=prctile(TKDc,options.subtract_prctile,1);
					tckd=tckd-plist;
					tckd(tckd<0)=0;
					CD=squeeze(mean(tckd,[1,3]));
				end
				score=sum(CD.*exp(1i*2*pi*[1:8]./8),2,"omitnan");	% cell x 1
				if options.divide_by_sum
					score=abs(score)./sum(CD,2);
				else
					score=abs(score);
				end
				[~,maxdir]=max(CD,[],2);
			end
			function pval=assess_ANOVA_pval(tckd)
				ckd=squeeze(mean(tckd,1));
				kdc=permute(ckd,[2,3,1]);
				pval=[];
				for i=1:size(kdc,3)
					pval(i,1)=anova1(kdc(:,:,i),[],'off');
				end
			end
		end
		function assess_direction_selectivity_use_cache_if_possible(obj,options)
			arguments
				obj
				options.use_margin_4_baseline (1,1) logical =false
				options.subtract_prctile (1,1)=[]
			end
			options=namedargs2cell(options);
			if isempty(obj.direction_selectivity)
				obj.assess_direction_selectivity(options{:});
			else
				if obj.plane~=obj.direction_selectivity.plane
					obj.assess_direction_selectivity(options{:});
				end
			end
		end
		%%
		function [cell_log,ROI]=plot_direction_selective_cell_density(obj,maskoptions,plotoptions,options,direction_selectivity)
			arguments
				obj
				maskoptions.mask_radius_in_um
				maskoptions.mask_shape (1,1) string {mustBeMember(maskoptions.mask_shape, ["Gaussian","Ellipse"])}
				maskoptions.padding_in_mm=0

				plotoptions.overlayColormap {mustBeNumeric} = jet(256)
				plotoptions.alphaVal (1,1) double {mustBeGreaterThanOrEqual(plotoptions.alphaVal,0), ...
					mustBeLessThanOrEqual(plotoptions.alphaVal,1)} = 1
				plotoptions.mask_thresh {mustBeNumeric} = []
				plotoptions.Title (1,:) char = ''
				plotoptions.imagescmax (1,1) double = 1
				plotoptions.Plot_mask4conv (1,1) logical=true

				options.selectivity_score_threshold=[]
				options.skewness_threshold=[]
				options.ANOVA_pval_thresh=[]
				options.show_plot =true
				options.plotHalf (1,1) logical =false

				direction_selectivity.pre_margin_in_sec=0
				direction_selectivity.post_margin_in_sec=0
				direction_selectivity.divide_by_sum (1,1) logical =false
				direction_selectivity.use_margin_4_baseline (1,1) logical =false
				direction_selectivity.subtract_prctile =[]
			end
			direction_selectivity=namedargs2cell(direction_selectivity);
			obj.assess_direction_selectivity(direction_selectivity{:});
			% if isempty(obj.direction_selectivity)
			% 	obj.assess_direction_selectivity(use_margin_4_baseline=options.use_margin_4_baseline);
			% else
			% 	if obj.plane~=obj.direction_selectivity.plane
				% 	obj.assess_direction_selectivity(use_margin_4_baseline=options.use_margin_4_baseline);
			% 	end
			% end
			ds=obj.direction_selectivity;
			score=ds.score;
			direction=ds.direction;
			if options.plotHalf
				ROI=obj.HalfROI;
			else
				ROI=obj.regROI;
			end
			skew_log=assess_high_skewness(obj.Skew,options.skewness_threshold);
			pval_log=assess_ANOVA_pval(ds.pval,options.ANOVA_pval_thresh);
			score_log=score>=options.selectivity_score_threshold;
			
			ROI=ROI(skew_log,:);
			pval_log=pval_log(skew_log);
			score_log=score_log(skew_log);

			cell_log=pval_log(:).*score_log(:);

			maskoptions=namedargs2cell(maskoptions);
			mask=CSP.generate_kernel_4_plotCellDensity(maskoptions{:});

			plotoptions=namedargs2cell(plotoptions);
			if options.show_plot
				if options.plotHalf
					plotCellDensity(ROI,cell_log,mask,plotoptions{:}, ...
						AtlasImage=obj.normal_HalfAtlas(),AtlasImage_mask=obj.make_HalfAtlas_mask());
				else
					plotCellDensity(ROI,cell_log,mask,plotoptions{:}, ...
						AtlasImage=obj.normal_AtlasImage(),AtlasImage_mask=obj.make_AtlasImage_mask());
				end
				title(sprintf("%s, # of cells: %d / %d (%.2g%%)",obj.PlaneName(obj.plane),sum(cell_log),numel(cell_log),100*sum(cell_log)/numel(cell_log)));
				drawnow;
			end
			
			function skewlog=assess_high_skewness(skew,skew_thresh)
				if isempty(skew_thresh)
					skewlog=true(numel(skew),1);
				else
					skewlog=skew>=skew_thresh;
				end
			end
			function pvallog=assess_ANOVA_pval(pval,pval_thresh)
				if isempty(pval_thresh)
					pvallog=true(numel(pval),1);
				else
					pvallog=pval<=pval_thresh;
				end
			end
		end

		%%
		function [all_cell_log,all_ROI]=plot_direction_selective_cell_density_all_planes_combined(obj,maskoptions,plotoptions,options,direction_selectivity)
			arguments
				obj
				maskoptions.mask_radius_in_um
				maskoptions.mask_shape (1,1) string {mustBeMember(maskoptions.mask_shape, ["Gaussian","Ellipse"])}
				maskoptions.padding_in_mm=0

				plotoptions.overlayColormap {mustBeNumeric} = jet(256)
				plotoptions.alphaVal (1,1) double {mustBeGreaterThanOrEqual(plotoptions.alphaVal,0), ...
					mustBeLessThanOrEqual(plotoptions.alphaVal,1)} = 1
				plotoptions.mask_thresh {mustBeNumeric} = []
				plotoptions.Title (1,:) char = ''
				plotoptions.imagescmax (1,1) double = 1
				plotoptions.Plot_mask4conv (1,1) logical=true

				options.selectivity_score_threshold=[]
				options.skewness_threshold=[]
				options.ANOVA_pval_thresh=[]
				options.show_plot =true
				options.plotHalf (1,1) logical =false

				direction_selectivity.pre_margin_in_sec=0
				direction_selectivity.post_margin_in_sec=0
				direction_selectivity.divide_by_sum (1,1) logical =false
				direction_selectivity.use_margin_4_baseline (1,1) logical =false
				direction_selectivity.subtract_prctile =[]
			end
			scorethresh=options.selectivity_score_threshold;
			skewthresh=options.skewness_threshold;
			ANOVAp=options.ANOVA_pval_thresh;
			options.show_plot=false;
			maskoptions=namedargs2cell(maskoptions);
			plotoptions=namedargs2cell(plotoptions);
			direction_selectivity=namedargs2cell(direction_selectivity);
			orgoptions=options;
			options=namedargs2cell(options);
			all_cell_log=[];
			all_ROI=[];
			for i=1:length(obj.csp_plane)
				csp=obj.csp_plane{i};
				[cell_log,ROI]=csp.plot_direction_selective_cell_density(maskoptions{:},plotoptions{:},options{:},direction_selectivity{:});
				all_cell_log=cat(1,all_cell_log,cell_log);
				all_ROI=cat(1,all_ROI,ROI);
				i
			end
			mask=CSP.generate_kernel_4_plotCellDensity(maskoptions{:});
			if orgoptions.plotHalf
				plotCellDensity(all_ROI,all_cell_log,mask,plotoptions{:},...
					AtlasImage=obj.normal_HalfAtlas(),AtlasImage_mask=obj.make_HalfAtlas_mask());
			else
				plotCellDensity(all_ROI,all_cell_log,mask,plotoptions{:},...
					AtlasImage=obj.normal_AtlasImage(),AtlasImage_mask=obj.make_AtlasImage_mask());
			end
			title(sprintf("Direction Selectivity of All Mice Combined. # of cells: %d / %d (%.2g%%), \nskew>=%.2g, ANOVA p<%.3g, Selectivity score>.%3g", ...
				sum(all_cell_log),numel(all_cell_log),100*sum(all_cell_log)/numel(all_cell_log),skewthresh,ANOVAp,scorethresh));
		end
	end
	%% direction selectivity
	methods
		function get_direction_selectivity(obj,options)
			arguments
				obj
				options.pre_margin_in_sec=0
				options.post_margin_in_sec=0
				options.divide_by_sum (1,1) logical =false
				options.use_margin_4_baseline (1,1) logical =false
				options.subtract_prctile =[]
				options.tckd=[]
				options.tpkd=[]
			end
			tckd=options.tckd;
			tpkd=options.tpkd;
			if isempty(options.tckd)
				tckd=obj.p8_tckd;
			end
			if isempty(options.tpkd)
				tpkd=obj.p8_tpkd;
			end



			if options.use_margin_4_baseline && ~isempty(options.subtract_prctile)
				warning("use_margin_4_baseline and subtract_prctile is incompatible")
			end
			obj.direction_selectivity=struct;
			tckd=obj.p8_tckd;
			pval=assess_ANOVA_pval(tckd);
			[score,maxdir,tckd]=assess_direction_selectivity(obj);
			obj.direction_selectivity.pval=pval;
			obj.direction_selectivity.score=score;
			obj.direction_selectivity.direction=maxdir;
			obj.direction_selectivity.plane=obj.plane;
			% obj.direction_selectivity.tckd=tckd;


			function [score,maxdir,tckd]=assess_direction_selectivity(obj)
				tckd=obj.p8_tckd;
				CD=squeeze(mean(tckd,[1,3]));
				csp=copy(obj);
				pre_mrgn=options.pre_margin_in_sec;
				post_mrgn=options.post_margin_in_sec;
				if pre_mrgn~=0 || post_mrgn~=0
					csp.loadData_with_trial_margin(pre_margin_in_sec=pre_mrgn,post_margin_in_sec=post_mrgn);
					csp.setCurrentData("p8");
					mtckd=csp.tck;
					pre_mrgn=pre_mrgn*obj.SamplingRate_Hz;
					post_mrgn=post_mrgn*obj.SamplingRate_Hz;
					base_tckd=mtckd([1:pre_mrgn,size(mtckd,1)-post_mrgn+1:size(mtckd,1)],:,:,:);
					base_tckd=mean(base_tckd,1);
					tckd=tckd-base_tckd;
					ckd=permute(mean(tckd,1),[2,3,4,1]);
					ckd(ckd<0)=0;
					CD=permute(mean(ckd,2),[1,3,2]);
				end
				if options.use_margin_4_baseline
					csp=copy(obj);
					mrgnsec=1;
					csp.loadData_with_trial_margin(pre_margin_in_sec=mrgnsec, ...
						post_margin_in_sec=mrgnsec);
					csp.setCurrentData("p8");
					mtckd=csp.tck;
					mrgn_duration=mrgnsec*obj.SamplingRate_Hz;
					base_tckd=mtckd([1:mrgn_duration,size(mtckd,1)-mrgn_duration:size(mtckd,1)], ...
						:,:,:);
					base_tckd=mean(base_tckd,1);
					tckd=tckd-base_tckd;
					tckd(tckd<0)=0;
					CD=squeeze(mean(tckd,[1,3]));
				end
				if ~isempty(options.subtract_prctile)
					TKDc=reshape(permute(tckd,[2,1,3,4]),size(tckd,2),[])';
					plist=prctile(TKDc,options.subtract_prctile,1);
					tckd=tckd-plist;
					tckd(tckd<0)=0;
					CD=squeeze(mean(tckd,[1,3]));
				end
				score=sum(CD.*exp(1i*2*pi*[1:8]./8),2,"omitnan");	% cell x 1
				if options.divide_by_sum
					score=abs(score)./sum(CD,2);
				else
					score=abs(score);
				end
				[~,maxdir]=max(CD,[],2);
			end
			function pval=assess_ANOVA_pval(tckd)
				ckd=squeeze(mean(tckd,1));
				kdc=permute(ckd,[2,3,1]);
				pval=[];
				for i=1:size(kdc,3)
					pval(i,1)=anova1(kdc(:,:,i),[],'off');
				end
			end
		end
	end
	%% generate mask 4 plotCellDensity
	methods (Static)
		function mask=generate_kernel_4_plotCellDensity(options)
			arguments
				options.mask_radius_in_um
				options.mask_shape (1,1) string {mustBeMember(options.mask_shape, ["Gaussian","Ellipse"])}
				options.padding_in_mm=0
			end
			mask_radius_in_um=options.mask_radius_in_um;
			mask_shape=options.mask_shape;
			padding_in_mm=options.padding_in_mm;
			one_cm=1.77;		% length of one centimeter in powerpoint
			width=18.89;		% length of image width in powerpoint
			height=15.02;		% length of image height in powerpoint
			pixwidth=911;
			pixheight=904;
			x_pix_per_mm=pixwidth*(one_cm/width);
			y_pix_per_mm=pixheight*(one_cm/height);
			mask=generate_mask4conv(x_pix_per_mm,y_pix_per_mm,mask_radius_in_um,mask_shape,padding_in_mm);
			obj.kernel_4_plotCellDensity=mask;
		end
	end
	methods (Static)
		function plotScaleBar(options)
			arguments
				options.BarLength_mm =1
				options.WidthPixel=10
				options.OriginPixel=[1,1]
				options.Color="black"
				options.FontSize=20
				options.ScaleXYequal (1,1) logical =true
				options.axis_display (1,1) logical =false
				options.axis_off=true
			end
			one_cm=1.77;		% length of one centimeter in powerpoint
			width=18.89;		% length of image width in powerpoint
			height=15.02;		% length of image height in powerpoint
			pixwidth=911;
			pixheight=904;
			x_pix_per_mm=pixwidth*(one_cm/width);
			y_pix_per_mm=pixheight*(one_cm/height);

			x_length=x_pix_per_mm*options.BarLength_mm;
			y_length=y_pix_per_mm*options.BarLength_mm;
			
			orgx=options.OriginPixel(1);
			orgy=options.OriginPixel(2);
			width=options.WidthPixel;

			% x bar
			x=[orgx,orgx,x_length+orgx,x_length+orgx;
				orgx,orgx,orgx+width,orgx+width]';
			y=[orgy,orgy+width,orgy+width,orgy;
				orgy,orgy+y_length,orgy+y_length,orgy]';
			patch(x,y,options.Color,EdgeColor="none");

			text(orgx+width,orgy+width,sprintf("%d mm",options.BarLength_mm),FontSize=options.FontSize,VerticalAlignment="top",HorizontalAlignment="left", ...
				Color=options.Color)

			if options.ScaleXYequal
				daspect([x_length/y_length,1,1]);
				% axis off
			end
			if ~options.axis_display
				if options.axis_off
					axis off
				end
			end
		end

	end
	%% generate density plot
	methods
		function plotCellDensity_all_plane(obj,options)
			arguments
				obj
				options.func_handle
			end
		end
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Ridge regression with time shift
	properties
		TimeShiftInSec_vec
	end
	methods
		function TimeShifted_Initialize_timeshift(obj,options)
			arguments
				obj
				options.TimeShiftInSec_vec
			end
			obj.TimeShiftInSec_vec=options.TimeShiftInSec_vec;
		end
		function TimeShifted_SetUpAndStartRidgeTraining(obj,options)
			% execute this after obj.loadData and obj.TimeShifted_Initialize_timeshift
			arguments
				obj
				options.ExecutionSequence_From_setCurrentData_to_Ridge_startTraining
			end
			InputList=["csp","plane","timeshift"];
			funcs=options.ExecutionSequence_From_setCurrentData_to_Ridge_startTraining;
			for i=1:length(obj.TimeShiftInSec_vec)
				tic
				obj.makeTimeShiftedData(TimeShiftInSec=obj.TimeShiftInSec_vec(i));
				executeFuncSequentially(obj,obj.plane,i,funcs);
				fprintf("%d / %d complete",i,length(obj.TimeShiftInSec_vec))
				toc
			end
			function executeFuncSequentially(csp,plane,timeshift,funcs)
				inputs={csp,plane,timeshift};
				for j=1:length(funcs)
					f=funcs{j};
					[nInputs,inputNames]=getFuncHandleInputs(f);
					inputind=arrayfun(@(a)find(InputList==a),inputNames);
					f(inputs{inputind});
				end
			end
		end
	end
%% Ridge regression
% Example:
% %% ridge rigression of all data
% groupnum=4;
% phibincount=5;
% DLC_bin_count=6;
% DLC_edge_prc_thresh=2;	% close to 0 %
% DLC_triError_thresh=3;	% threshold for triagulation error (mm)
% sigma=1;
% ridgeparams=10.^[-10:10];
% 
% ridgeWeights=[];
% ridgeResults=[];
% 
% 
% allstart=tic;
% displayer = Progressdisp(length(info.NoAltplane.nameMouse_p), barWidth=40, BarOrStr="bar");
% for plane=1:length(info.NoAltplane.nameMouse_p)
% 	Data=CSP(info=info);
% 	Data.loadData(plane);
% 	Data.setCurrentData("p8");
% 	Data.GaussHotEncode_SPparams(sigma,{makeThetaBins,makePhiBins(phibincount)});
% 	Data.Generate_edge4gausDLC(DLC_triError_thresh,DLC_bin_count,DLC_edge_prc_thresh,"invDLC");
% 	Data.GaussHotEncode_DLC(sigma,["tXYZbp","tDXYZbp","tD2XYZbp","tANGbp","tDANGbp","tD2ANGbp"]);
% 	Data.combinePredictors("Both");
% 	Data.Ridge_makeCrossValidationGroups(groupnum);
% 	ridgeResults{plane}=Data.Ridge_startTraining(ridgeparams);
% 	displayer.Display_progress(plane);
% end
% toc(allstart);




% Gaussian-hot encoding
	properties
		tpk_bidirectional
		EncodeBidirectionalPhi
		SP_GaussianHotEncoded_cellVector
		SP_GaussianHotEncoded_Discretizer
		SP_GaussianHotEncoded_Undiscretizer
		SP_GaussianHotEncoded_edges
		SP_GaussianHotEncoded_bincount
		SP_GaussianHotEncoded_sigma
		SP_GaussianHotEncoded_params2encode
		SP_GaussianHotEncoded_params2encode_cell

		SP_GaussianHotEncoded_names

		DLC_GaussianHotEncoded_cellVector
		DLC_GaussianHotEncoded_edges_cellVector
		DLC_GaussianHotEncoded_sigma
		DLC_GaussianHotEncoded_params2encode

		all_GaussianHot
		all_GaussianHot_option
		all_GaussianHot_names


		SP_edges4GaussianHotEncoding

		gausHot_SP_sigma	% sigma of gaussian for SP
		edge4gausSP
		gausHot_DLC_sigma	% sigma of gaussian for DLC
		DLC_triError_thresh4EdgeCalculation
		DLC_edge_prc_thresh4EdgeCalculation
		edge4gausDLC
		DLC_option
		gausHot_SP
		gausHot_DLC
		cat_gausHot_SP
		cat_gausHot_DLC
		all_gausHot
	end


	methods
		function SP_GaussianHotEncode(obj,sigma,params2discretize_stringArray4jointDiscretization,options)
			arguments
				obj
				sigma
			end
			arguments (Repeating)
				params2discretize_stringArray4jointDiscretization (1,:) string {mustBeMember(params2discretize_stringArray4jointDiscretization, ...
					["theta","phi","tau","X","Y","Z"])}
			end
			arguments
				options.Output (1,1) string {mustBeMember(options.Output, ["CellArray","NumericArray_concatenateOutputIn2ndDimension"])}="CellArray"
				options.EncodeBidirectionalPhi logical =true
			end
			if options.EncodeBidirectionalPhi
				obj.EncodeBidirectionalPhi=true;
				[newtpk,newedges]=make_tpk_and_edges_alongwith_bidirectional_phi();
				tpk=newtpk;
				edges=newedges;
				obj.tpk_bidirectional=tpk;
			else
				obj.EncodeBidirectionalPhi=false;
				tpk=obj.tpk;
				edges=obj.SP_edges_cellVector;
			end
			options=rmfield(options,"EncodeBidirectionalPhi");
			options=namedargs2cell(options);
			[GaussHotEncoded_Data,Func]=GaussianHotEncode_SpecifedSPparams...
				(tpk,edges,sigma,params2discretize_stringArray4jointDiscretization{:}, options{:});
			Func=CellWithStruct2StructWithCell(Func);

			obj.SP_GaussianHotEncoded_cellVector=GaussHotEncoded_Data;
			% make names
			obj.SP_GaussianHotEncoded_names=[];
			names={};
			for i=1:length(params2discretize_stringArray4jointDiscretization)
				temp=params2discretize_stringArray4jointDiscretization{i};
				temp=strjoin(temp,"_");
				temp=repmat(temp,1,size(obj.SP_GaussianHotEncoded_cellVector{i},2));
				obj.SP_GaussianHotEncoded_names{i}=temp;
			end
			if assessNaN(GaussHotEncoded_Data)
				warning(sprintf("There are NaN in GaussianHotEncoded_Data of plane %d",obj.plane))
			end
			obj.SP_GaussianHotEncoded_Discretizer=Func.discretize;
			obj.SP_GaussianHotEncoded_Undiscretizer=Func.undiscretize;
			obj.SP_GaussianHotEncoded_edges=Func.edges;
			obj.SP_GaussianHotEncoded_bincount=Func.bincount;
			obj.SP_GaussianHotEncoded_sigma=sigma;
			obj.SP_GaussianHotEncoded_params2encode=[params2discretize_stringArray4jointDiscretization{:}];
			obj.SP_GaussianHotEncoded_params2encode_cell=params2discretize_stringArray4jointDiscretization;
			function [newtpk,newedges]=make_tpk_and_edges_alongwith_bidirectional_phi()
				IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam("phi");
				currdata=obj.tpk(IND{:});
				orgdata=currdata;
				diffdata=diff(currdata,1,1);
				IND=repmat({':'},1,ndims(obj.tpk));IND{1}=size(diffdata,1);
				altIND=IND;altIND{1}=size(diffdata,1)+1;
				diffdata(altIND{:})=diffdata(IND{:});		% make the size same as currdata
				currdata(diffdata==0)=NaN;
				currdata(diffdata>0)=currdata(diffdata>0);
				currdata(diffdata<0)=-currdata(diffdata<0);	% ascending --> positive, descending --> negative
				currdata=fillmissing(currdata,"previous",1);
				currdata=fillmissing(currdata,"next",1);
				currdata=replaceAllNaN_in1stdimension_toOriginal(currdata,orgdata);
				
				phiedges=obj.SP_edges_cellVector{obj.ind_SPparam("phi")};
				nphiedges=-phiedges;
				phiedges=[nphiedges(end:-1:2);phiedges];

				% assign new variables
				newtpk=obj.tpk;IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam("phi");
				newtpk(IND{:})=currdata;
				newedges=obj.SP_edges_cellVector;
				newedges{obj.ind_SPparam("phi")}=phiedges;
			end
			function outputdata=replaceAllNaN_in1stdimension_toOriginal(nanned_data,orginal_data)
				nancheck=all(isnan(nanned_data),1);
				nanlog=repmat(nancheck,size(nanned_data,1),1);
				nanned_data(nanlog)=orginal_data(nanlog);
				outputdata=nanned_data;
			end
		end
		%% gaussian hot encode using Discretized SP
		function SP_GaussianHotEncode_Discretized(obj,sigma)
			arguments
				obj
				sigma
			end
			encoded_data={};
			for i=1:numel(obj.SP_Discretized_cellVector)
				data=obj.SP_Discretized_cellVector{i};
				bins=obj.SP_Discretized_bincount{i};
				encoded_data{i}=GaussianHotEncoder_from_Discretized(data,1,2,sigma,bins);
			end
			obj.SP_GaussianHotEncoded_cellVector=encoded_data;
			% make names
			obj.SP_GaussianHotEncoded_names=[];
			names={};
			for i=1:length(obj.SP_Discretized_params2encode_cell)
				temp=obj.SP_Discretized_params2encode_cell{i};
				temp=strjoin(temp,"_");
				temp=repmat(temp,1,size(obj.SP_GaussianHotEncoded_cellVector{i},2));
				obj.SP_GaussianHotEncoded_names{i}=temp;
			end

			obj.SP_GaussianHotEncoded_Discretizer=obj.SP_Discretized_Discretizer;
			obj.SP_GaussianHotEncoded_Undiscretizer=obj.SP_Discretized_Undiscretizer;
			obj.SP_GaussianHotEncoded_edges=obj.SP_Discretized_edges;
			obj.SP_GaussianHotEncoded_bincount=obj.SP_Discretized_bincount;
			obj.SP_GaussianHotEncoded_sigma=sigma;
			obj.SP_GaussianHotEncoded_params2encode=obj.SP_Discretized_params2encode;
			obj.SP_GaussianHotEncoded_params2encode_cell=obj.SP_Discretized_params2encode_cell;
		end
		%%
		function plot_SP_GaussianHotEncoded(obj,param2plot)
			f=@(a)BlockCat(reshape(a,size(a,1),size(a,2),[]),3,1);
			ax1=nexttile;
			currmotions=obj.SP_GaussianHotEncoded_params2encode_cell{param2plot};
			IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam(currmotions);
			data=obj.tpk(IND{:});
			data=f(data);
			plot(data);

			ax2=nexttile;
			data=obj.SP_GaussianHotEncoded_cellVector{param2plot};
			data=f(data);
			imagesc(data')
			linkaxes([ax1,ax2],'x');
		end
		% function plot_SP_GaussianHotEncoded_alt(obj,param2plot)
		% 	f=@(a)BlockCat(reshape(a,size(a,1),size(a,2),[]),3,1);
		% 	ax1=nexttile;
		% 	currmotions=obj.SP_GaussianHotEncoded_params2encode_cell{param2plot};
		% 	IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam(currmotions);
		% 	data=obj.tpk(IND{:});
		% 	data=f(data);
		% 	plot(data);
		% 
		% 	ax2=nexttile;
		% 	data=obj.SP_GaussianHotEncoded_cellVector{param2plot};
		% 	data=f(data);
		% 	imagesc(data')
		% 	linkaxes([ax1,ax2],'x');
		% end
		
		%%
		function plot_SP_GaussianHotEncoded_new(obj,param2plot)
			f=@(a)BlockCat(reshape(a,size(a,1),size(a,2),[]),3,1);
			ax1=nexttile;
			currmotions=obj.SP_GaussianHotEncoded_params2encode_cell{param2plot};
			IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam(currmotions);
			data=obj.tpk(IND{:});
			data=f(data);
			plot(data);

			ax2=nexttile;
			data=obj.SP_GaussianHotEncoded_cellVector{param2plot};
			data=f(data);
			imagesc(data')
			linkaxes([ax1,ax2],'x');
		end
	end
	properties
		DLC_GaussianHotEncoded_names
	end
	methods
		function DLC_GaussianHotEncode(obj,sigma,params2discretize)
			arguments
				obj
				sigma
			end
			arguments (Repeating)
				params2discretize (1,1) string {mustBeMember(params2discretize,["tXYZbp","tDXYZbp","tD2XYZbp", ...
					"tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp", ...
					"tDXYZbp_magnitude","tDDOTbp_magnitude","tDANGbp_magnitude", ...
					"tD2XYZbp_magnitude","tD2DOTbp_magnitude","tD2ANGbp_magnitude"])}
			end
			XYZ=["tXYZbp","tDXYZbp","tD2XYZbp"];
			Ang=["tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"];

			obj.DLC_GaussianHotEncoded_cellVector=cell(1,length(params2discretize));
			obj.DLC_GaussianHotEncoded_names=cell(1,length(params2discretize));
			switch obj.DLC_inv_or_normal
				case "normalPosture"
					DLCdata=obj.DLC_convert(obj.DLC,"AllNumeric");
				case "PostureInversedBySP"
					DLCdata=obj.DLC_convert(obj.invDLC,"AllNumeric");
			end
			for i=1:length(params2discretize)
				currstr=params2discretize{i};
				name_of_current_category=[];
				if isfield(DLCdata,currstr)
					currdata=DLCdata.(currstr);
					curredges=obj.DLC_edges_type_cellVector.(currstr);
					if any(containsSubstring(XYZ,currstr))
						temp=arrayfun(@(a)currstr+"_"+a+["_x","_y","_z"],obj.rec3d.bp_ind,'UniformOutput',false);
						name_of_current_category=cat(2,temp{:});...1 x 42 string array
					elseif any(containsSubstring(Ang,currstr))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif isfield(obj.DLC_speed,currstr)
					currdata=obj.DLC_speed.(currstr);
					curredges=obj.DLC_edges_speed_cellVector.(currstr);
					if any(containsSubstring(currstr,"XYZ"))
						name_of_current_category=currstr+"_"+obj.rec3d.bp_ind;...1 x 14 string array
					elseif any(containsSubstring(currstr,"ANG"))||any(containsSubstring(currstr,"DOT"))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif isfield(obj.DLC_accel,currstr)
					currdata=obj.DLC_accel.(currstr);
					curredges=obj.DLC_edges_accel_cellVector.(currstr);
					if any(containsSubstring(currstr,"XYZ"))
						name_of_current_category=currstr+"_"+obj.rec3d.bp_ind;...1 x 14 string array
					elseif any(containsSubstring(currstr,"ANG"))||any(containsSubstring(currstr,"DOT"))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				else
					error("No match")
				end
				obj.DLC_GaussianHotEncoded_cellVector{i}=GaussianHotEncoder(currdata,1,2,sigma,curredges);
				obj.DLC_GaussianHotEncoded_edges_cellVector{i}=curredges;
				name_of_current_category=cellstr(name_of_current_category);
				name_of_current_category=cellfun(@(a,b)repmat(string(a),1,numel(b)-1),name_of_current_category,curredges,UniformOutput=false);
				obj.DLC_GaussianHotEncoded_names{i}=name_of_current_category;
			end

			% obj.DLC_GaussianHotEncoded_cellVector=cat(2,obj.DLC_GaussianHotEncoded_cellVector{:});
			% obj.DLC_GaussianHotEncoded_edges_cellVector=cat(2,obj.DLC_GaussianHotEncoded_edges_cellVector{:});
			obj.DLC_GaussianHotEncoded_sigma=sigma;
			obj.DLC_GaussianHotEncoded_params2encode=[params2discretize{:}];
		end
		%%
		function plot_DLC_GaussianHotEncoded(obj,params2plot,ind)
			arguments
				obj
				params2plot (1,1) string {mustBeMember(params2plot,["tXYZbp","tDXYZbp","tD2XYZbp","tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp", ...
					"tDXYZbp_magnitude","tDDOTbp_magnitude","tDANGbp_magnitude", ...
					"tD2XYZbp_magnitude","tD2DOTbp_magnitude","tD2ANGbp_magnitude"])}
				ind
			end
			currstrIND=find(obj.DLC_GaussianHotEncoded_params2encode==params2plot);
			switch obj.DLC_inv_or_normal
				case "normalPosture"
					DLCdata=obj.DLC_convert(obj.DLC,"AllNumeric");
				case "PostureInversedBySP"
					DLCdata=obj.DLC_convert(obj.invDLC,"AllNumeric");
			end
			try
				tpk=DLCdata.(params2plot);
			catch
				try
					tpk=obj.DLC_speed.(params2plot);
				catch
					tpk=obj.DLC_accel.(params2plot);
				end
			end
			tpk_GaussHot=obj.DLC_GaussianHotEncoded_cellVector{currstrIND};
			edges=obj.DLC_GaussianHotEncoded_edges_cellVector{currstrIND};
			plotGaussianHotEncodedData(tpk,tpk_GaussHot,edges, ind)
		end
		%%
		function plot_DLC_GaussianHotEncoded_alt(obj,param2plot,ntl)
			arguments
				obj
				param2plot (1,1) string
				ntl=[]
			end

			for i=1:length(obj.DLC_Discretized_names)
				ind=obj.DLC_Discretized_names{i}==param2plot;
				if sum(ind)>0
					discretized=obj.DLC_GaussianHotEncoded_cellVector{i}{ind};
					undiscretized=obj.DLC_Discretized_Undiscretized{i}{ind};
					org_data=obj.DLC_Discretized_orgData{i}{ind};
					edges=obj.DLC_Discretized_edges_cellVector{i}{ind};
				end
			end
			f=@(a)multiBlockCat(a,[3:ndims(a)],1);
			plotHotEncoded(f(discretized),f(org_data),edges,ntl);
		end
		%%
		function DLC_GaussianHotEncode_Discretized(obj,sigma)
			arguments
				obj
				sigma
			end
			encoded_data={};
			names=obj.DLC_Discretized_names;
			bins=obj.DLC_Discretized_bincount;
			for i=1:numel(obj.DLC_Discretized_cellVector)
				data=obj.DLC_Discretized_cellVector{i};
				encoded_data{i}=cellfun(@(array)GaussianHotEncoder_from_Discretized(array,1,2,sigma,bins),data,UniformOutput=false);
				names{i}=cellstr(names{i});
				names{i}=cellfun(@(a)repmat(string(a),1,bins),names{i},UniformOutput=false);
				% encoded_data{i}=GaussianHotEncoder_from_Discretized(data,1,2,sigma,bins);
			end
			obj.DLC_GaussianHotEncoded_cellVector=encoded_data;
			obj.DLC_GaussianHotEncoded_edges_cellVector=obj.DLC_Discretized_edges_cellVector;
			obj.DLC_GaussianHotEncoded_names=names;
			obj.DLC_GaussianHotEncoded_sigma=sigma;
			obj.DLC_GaussianHotEncoded_params2encode=obj.DLC_Discretized_params2encode;
		end
		
%%
		function combinePredictors(obj,option,UseValuesAsPredictors,options)
			arguments
				obj
				option (1,1) string {mustBeMember(option, ["Both","DLC only","SP only"])}
				UseValuesAsPredictors (1,1) logical =true
				options.removeSmilarParameters (1,1) logical =false
			end
			SPdata=obj.SP_GaussianHotEncoded_cellVector;
			SPdata=cell2mat(SPdata);
			DLCdata=obj.DLC_GaussianHotEncoded_cellVector;
			DLCdata=cell2mat([DLCdata{:}]);
			% categorize names
			SPdata_names=obj.SP_GaussianHotEncoded_names;
			SPdata_names=cat(2,SPdata_names{:});
			DLCdata_names=obj.DLC_GaussianHotEncoded_names;
			temp=[DLCdata_names{:}];
			DLCdata_names=cat(2,temp{:});
			if size(SPdata,2)~=size(SPdata_names,2);error("The size of array and names do not match");end
			if size(DLCdata,2)~=size(DLCdata_names,2);error("The size of array and names do not match");end

			XYZ=["tXYZbp","tDXYZbp","tD2XYZbp"];
			Ang=["tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"];


			if UseValuesAsPredictors
				% SP data
				currstr=obj.SP_GaussianHotEncoded_params2encode;
				currstr=unique(currstr);
				IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam(currstr);
				SPdata=cat(2,SPdata,obj.tpk(IND{:}));
				% concatenate names
				currname=currstr;
				SPdata_names=cat(2,SPdata_names,currname);
				if size(SPdata,2)~=size(SPdata_names,2);error("The size of array and names do not match");end
				
				% DLC data
				switch obj.DLC_inv_or_normal
					case "normalPosture"
						currdata=obj.DLC_convert(obj.DLC,"AllNumeric");
					case "PostureInversedBySP"
						currdata=obj.DLC_convert(obj.invDLC,"AllNumeric");
				end

				DLC_normal=currdata;
				currstr=obj.DLC_GaussianHotEncoded_params2encode;
				for k=currstr(:)'
					name_of_current_category=[];
					if isfield(DLC_normal,k)
						DLCdata=cat(2,DLCdata,DLC_normal.(k));
						if any(containsSubstring(XYZ,k))
							temp=arrayfun(@(a)k+"_"+a+["_x","_y","_z"],obj.rec3d.bp_ind,'UniformOutput',false);
							name_of_current_category=cat(2,temp{:});...1 x 42 string array
						elseif any(containsSubstring(Ang,k))
							name_of_current_category=k+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
						else
							error("No match")
						end
					elseif isfield(obj.DLC_speed,k)
						temp=obj.DLC_speed;
						DLCdata=cat(2,DLCdata,temp.(k));
						if any(containsSubstring(k,"XYZ"))
							name_of_current_category=k+"_"+obj.rec3d.bp_ind;...1 x 14 string array
						elseif any(containsSubstring(k,"ANG"))||any(containsSubstring(k,"DOT"))
							name_of_current_category=k+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
						else
							error("No match")
						end
					elseif isfield(obj.DLC_accel,k)
						temp=obj.DLC_accel;
						DLCdata=cat(2,DLCdata,temp.(k));
						if any(containsSubstring(k,"XYZ"))
							name_of_current_category=k+"_"+obj.rec3d.bp_ind;...1 x 14 string array
						elseif any(containsSubstring(k,"ANG"))||any(containsSubstring(k,"DOT"))
							name_of_current_category=k+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
						else
							error("No match")
						end
					else
						error("No match")
					end
					name_of_current_category=name_of_current_category+"_raw";
					DLCdata_names=cat(2,DLCdata_names,name_of_current_category);
					if size(DLCdata,2)~=size(DLCdata_names,2);error("The size of array and names do not match");end
					% try
					% 	DLCdata=cat(2,DLCdata,currdata.(k));
					% catch
					% 	try
						% 	temp=obj.DLC_speed;
						% 	DLCdata=cat(2,DLCdata,temp.(k));
					% 	catch
						% 	temp=obj.DLC_accel;
						% 	DLCdata=cat(2,DLCdata,temp.(k));
					% 	end
					% end
				end
			end
			if size(SPdata,2)~=size(SPdata_names,2);error("The size of array and names do not match");end
			if size(DLCdata,2)~=size(DLCdata_names,2);error("The size of array and names do not match");end
			switch option
				case "Both"
					all_Pred=[SPdata,DLCdata];
					all_names=[SPdata_names,DLCdata_names];
				case "DLC only"
					all_Pred=[DLCdata];
					all_names=[DLCdata_names];
				case "SP only"
					all_Pred=[SPdata];
					all_names=[SPdata_names];
			end

			% remove predictors that are constant	%%%%%%%% not compatible
			% currently, names need to be adjusted if you wanna remove
			% similar parameters
			if options.removeSmilarParameters
				a=permute(all_Pred,[2,1,3:ndims(all_Pred)]);
				orgsize=size(a);
				X=reshape(a,size(a,1),[])';
				[X_new, keep_idx] = remove_duplicate_columns_from_matrix(X);
				rangelog=removeConstantColumns(X);
				keep_ind=intersect(keep_idx,rangelog);
				X_new=X(:,keep_ind);
				orgsize(1)=size(X_new,2);
				X_new=reshape(X_new',orgsize);
				X_new=permute(X_new,[2,1,3:ndims(X_new)]);
				all_Pred=X_new;
			end

			% normalize predictors
			IND=repmat({':'},1,ndims(all_Pred));
			sp_mu=mean(all_Pred,[1,3:ndims(all_Pred)]);
			sp_sigma=std(all_Pred,0,[1,3:ndims(all_Pred)]);
			% for constant variables, just leave them as they are
			sp_sigma(sp_sigma==0)=1;
			obj.all_GaussianHot=(all_Pred-sp_mu)./sp_sigma;
			obj.all_GaussianHot_option=option;
			obj.all_GaussianHot_names=all_names;

			function rangelog=removeConstantColumns(data)
				datarange=range(data,1);
				rangelog=find(datarange>0);
			end
		end
	end

	% Ridge regression
	properties
		B				% calculated weight
		triErrorThresh4DataOmission_ridge
		shiftIn10Hz_vec
		ridgeparams
		ridgeResults_CellRidgeparams	% cell x ridge parameters
		ridgeResults_setting_plane_CellRidge {mustBeA(ridgeResults_setting_plane_CellRidge,"struct")}=struct	% cell x ridge parameters array concatenated in cell vector
		ridgeResultsPlaneIND_setting_plane_CellRidge  {mustBeA(ridgeResultsPlaneIND_setting_plane_CellRidge,"struct")}=struct
		ridgeResults_setting_plane_CellRidgeShift =struct
	end
	methods

		%%
		function Ridge_makeCrossValidationGroups(obj,options)
			arguments
				obj
				options.groupnum
				options.triErrorThresh4DataOmission=[]
			end
			obj.triErrorThresh4DataOmission_ridge=options.triErrorThresh4DataOmission;
			groupnum=options.groupnum;
			if any(isnan(obj.all_GaussianHot(:)))
				error("There are NaN values in gaussianhot encoded data")
			end
			if isempty(options.triErrorThresh4DataOmission)
				pred=obj.all_GaussianHot;
				resp=obj.tck;
			else
				pred=obj.all_GaussianHot;
				triError=obj.DLC.triError;
				trilog=triError<=options.triErrorThresh4DataOmission;	% log for accurate time points
				trilog=prod(trilog,2);
				btrilog=trilog+zeros(size(pred));	% broadcast to all the predictor parameters (2nd dim)
				btrilog=~logical(btrilog);				% log for inaccurate timepoints
				pred(btrilog)=NaN;
				resp=obj.tck;
				btrilog=trilog+zeros(size(resp));
				btrilog=~logical(btrilog);				% log for inaccurate timepoints
				resp(btrilog)=NaN;
			end
			% split into cross validation groups
			[Trainer, Tester] = splitPredResp(pred, resp, groupnum,3);
			f=@reshapeAndConcatenate3rdTo1stDim;	% concatenate all 4th and higher dimensions to 3rd, and concatenate 3rd to 1st dimension to obtain 2d matrix
			obj.Trainer=applyFunctionLayers(Trainer,f);
			obj.Tester=applyFunctionLayers(Tester,f);
			obj.groupnum=groupnum;

			if isempty(options.triErrorThresh4DataOmission)
			else
				f=@(a)RemoveNanDataPoints(a);
				obj.Trainer=applyFunctionLayers(obj.Trainer,f);
				obj.Tester=applyFunctionLayers(obj.Tester,f);
				Log1=SearchEmptyElementFromAnyStruct(obj.Trainer);
				Log2=SearchEmptyElementFromAnyStruct(obj.Tester);
				Log=Log1.*Log2;
				f=@(a)a(Log);
				% f=@removeEmptyElementFromCell;
				obj.Trainer=structfun(f,obj.Trainer,UniformOutput=false);
				obj.Tester=structfun(f,obj.Tester,UniformOutput=false);
			end
			function newdata=RemoveNanDataPoints(data)
				if ndims(data)>2
					error("Data must be a matrix");
				end
				if sum(abs(diff(isnan(data),1,2)),"all")>0
					error("Each row must be either all NaN or all non-NaN");
				end
				a=data(:,1);
				newdata=data(~isnan(a),:);
				if any(isnan(newdata(:)))
					error("NaN value remains inside matrix");
				end
			end
			function emptylog=SearchEmptyElementFromAnyStruct(data)
				celllog=[];
				flds1=string(fields(data));
				for i=flds1(:)'
					temp=~cellfun(@isempty,data.(i),UniformOutput=true);
					celllog=cat(1,celllog,temp);
					if any(temp(:)==0,"all")
						warning("There are empty cells");
					end
				end
				emptylog=prod(celllog,1);
			end
			function newdata=removeEmptyElementFromCell(celldata)
				if ~isvector(celldata)
					error("celldata must be a vector cell");
				end
				newdata=celldata(~cellfun(@isempty,celldata,UniformOutput=true));
			end
			function newdata=reshapeAndConcatenate3rdTo1stDim(data)
				data=reshape(data,size(data,1),size(data,2),[]);
				data=permute(data,[2,1,3]);		% cell x time x trials
				newdata=reshape(data,size(data,1),[])';	% (cell x time*trials)'
			end
		end
		%% make xvalidation groups for continuous data
		function Ridge_makeCrossValidationGroups_4continuous(obj,options)
			arguments
				obj
				options.groupnum
				options.trialInterval_in_sec=[]
				options.triErrorThresh4DataOmission=[]
			end
			obj.triErrorThresh4DataOmission_ridge=options.triErrorThresh4DataOmission;
			groupnum=options.groupnum;
			if any(isnan(obj.all_GaussianHot(:)))
				error("There are NaN values in gaussianhot encoded data")
			end
			if isempty(options.triErrorThresh4DataOmission)
				pred=obj.all_GaussianHot;
				resp=obj.tck;
			else
				pred=obj.all_GaussianHot;
				triError=obj.DLC.triError;
				trilog=triError<=options.triErrorThresh4DataOmission;	% log for accurate time points
				trilog=prod(trilog,2);
				btrilog=trilog+zeros(size(pred));	% broadcast to all the predictor parameters (2nd dim)
				btrilog=~logical(btrilog);				% log for inaccurate timepoints
				pred(btrilog)=NaN;
				resp=obj.tck;
				btrilog=trilog+zeros(size(resp));
				btrilog=~logical(btrilog);				% log for inaccurate timepoints
				resp(btrilog)=NaN;
			end
			% split data to 3 dimensions
			if size(pred,1)~=size(resp,1)
				error("Size of pred and resp must match");
			end
			pred=reshapeTo3rdDim(pred,options.groupnum);
			resp=reshapeTo3rdDim(resp,options.groupnum);

			% split into cross validation groups
			[Trainer, Tester] = splitPredResp(pred, resp, groupnum,3);
			f=@reshapeAndConcatenate3rdTo1stDim;	% concatenate all 4th and higher dimensions to 3rd, and concatenate 3rd to 1st dimension to obtain 2d matrix
			obj.Trainer=applyFunctionLayers(Trainer,f);
			obj.Tester=applyFunctionLayers(Tester,f);
			obj.groupnum=groupnum;

			if isempty(options.triErrorThresh4DataOmission)
			else
				f=@(a)RemoveNanDataPoints(a);
				obj.Trainer=applyFunctionLayers(obj.Trainer,f);
				obj.Tester=applyFunctionLayers(obj.Tester,f);
				Log1=SearchEmptyElementFromAnyStruct(obj.Trainer);
				Log2=SearchEmptyElementFromAnyStruct(obj.Tester);
				Log=Log1.*Log2;
				f=@(a)a(Log);
				% f=@removeEmptyElementFromCell;
				obj.Trainer=structfun(f,obj.Trainer,UniformOutput=false);
				obj.Tester=structfun(f,obj.Tester,UniformOutput=false);
			end
			function newmat=reshapeTo3rdDim(mat,groupnum)
				timelength=floor(size(mat,1)/groupnum);
				mat=mat(1:timelength*groupnum,:);
				newmat=permute(reshape(mat',size(mat,2),timelength,groupnum),[2,1,3]);
			end
			function newdata=RemoveNanDataPoints(data)
				if ndims(data)>2
					error("Data must be a matrix");
				end
				if sum(abs(diff(isnan(data),1,2)),"all")>0
					error("Each row must be either all NaN or all non-NaN");
				end
				a=data(:,1);
				newdata=data(~isnan(a),:);
				if any(isnan(newdata(:)))
					error("NaN value remains inside matrix");
				end
			end
			function emptylog=SearchEmptyElementFromAnyStruct(data)
				celllog=[];
				flds1=string(fields(data));
				for i=flds1(:)'
					temp=~cellfun(@isempty,data.(i),UniformOutput=true);
					celllog=cat(1,celllog,temp);
					if any(temp(:)==0,"all")
						warning("There are empty cells");
					end
				end
				emptylog=prod(celllog,1);
			end
			function newdata=removeEmptyElementFromCell(celldata)
				if ~isvector(celldata)
					error("celldata must be a vector cell");
				end
				newdata=celldata(~cellfun(@isempty,celldata,UniformOutput=true));
			end
			function newdata=reshapeAndConcatenate3rdTo1stDim(data)
				data=reshape(data,size(data,1),size(data,2),[]);
				data=permute(data,[2,1,3]);		% cell x time x trials
				newdata=reshape(data,size(data,1),[])';	% (cell x time*trials)'
			end
		end
		%%
		function Ridge_makeCrossValidationGroups_ShiftTimePointsOfCellData(obj,options)
			arguments
				obj
				options.groupnum
				options.triErrorThresh4DataOmission=[]
				options.shiftIn10Hz_vec=0;
			end
			obj.shiftIn10Hz_vec=options.shiftIn10Hz_vec;
			obj.groupnum=options.groupnum;
			options=namedargs2cell(rmfield(options,"shiftIn10Hz_vec"));
			orgtck=obj.tck;
			obj.Trainer=[];
			obj.Tester=[];
			for i=1:length(obj.shiftIn10Hz_vec)
				csp=copy(obj);
				shift=obj.shiftIn10Hz_vec(i);
				[csp.all_GaussianHot,csp.tck]=makeShiftIND(obj.all_GaussianHot,obj.tck,shift);	% replace data with the shifted ones
				[csp.DLC.triError,~]=makeShiftIND(csp.DLC.triError,obj.tck,shift);				% replace data with the shifted ones
				csp.Ridge_makeCrossValidationGroups(options{:});
				obj.Trainer{i}=csp.Trainer;
				obj.Tester{i}=csp.Tester;
			end
			function [shiftedPred,shiftedResp]=makeShiftIND(pred,resp,shift)
				p=size(pred,1);
				r=size(resp,1);	% ==p
				[pind,rind]=ShiftAndGetOverlappingVec(1:p,1:r,shift);
				IND=createIND(pred,1,pind);
				shiftedPred=pred(IND{:});
				IND=createIND(resp,1,rind);
				shiftedResp=resp(IND{:});
			end
		end

		function Ridge_makeCrossValidationGroups_TimeShift(obj,options)
			arguments
				obj
				options.groupnum
				options.triErrorThresh4DataOmission=[]
				options.shiftIn10Hz_vec=0;
			end
			options
		end
	end
	%% 
	properties
		nameOfCurrentSettingRidge
		Meanscore_CellRidgeparams
		Stdscore_CellRidgeparams
		Weights_ridge
		Resp_group
		Rec_group
	end
	methods
		%%
		function ridgeResults=Ridge_startTraining(obj,ridgeparams,nameOfCurrentSetting,options)
			arguments
				obj
				ridgeparams
				nameOfCurrentSetting (1,1) string
				options.MakeDataSingle=true
				options.DisplayGPUMemory=false
				options.Directory2saveData=[]
				options.DisplayTimeTook2Save (1,1) logical =false
				options.SaveDemoCellActivity (1,1) logical =false
				options.SaveWeights	(1,1) logical =false
				options.saveCellActivity (1,1) logical =false
			end
			obj.nameOfCurrentSettingRidge=nameOfCurrentSetting;

			Trainer=obj.Trainer;
			Tester=obj.Tester;
			currResults=[];
			obj.B=[];
			obj.Rec=[];
			B_group=[];
			y_group=[];
			y_recon_group=[];
			for i=1:obj.groupnum
				if options.DisplayGPUMemory
					fprintf("Available GPU Memory: %d\n", gpuDevice().AvailableMemory);
				end
				 [currResults{i,1},B_group{i},y_group{i},y_recon_group{i}] = processGroup(Trainer, Tester, ridgeparams, options.MakeDataSingle, i);
			end
			currResults=cell2mat(currResults);
			ridgeResults=DeleteSingleDim(mean(currResults,1),1);	% cell x ridge parameters
			obj.ridgeResults_CellRidgeparams=ridgeResults;
			obj.ridgeparams=ridgeparams;

			if ~isfield(obj.ridgeResults_setting_plane_CellRidge,nameOfCurrentSetting)
				obj.ridgeResults_setting_plane_CellRidge.(nameOfCurrentSetting)=cell(0);
				obj.ridgeResultsPlaneIND_setting_plane_CellRidge.(nameOfCurrentSetting)=cell(0);
			end
			obj.ridgeResults_setting_plane_CellRidge.(nameOfCurrentSetting){end+1}=ridgeResults;
			obj.ridgeResultsPlaneIND_setting_plane_CellRidge.(nameOfCurrentSetting){end+1}=obj.plane;

			[Meanscore,Stdscore]=calculateRegressionPerformance(y_group,y_recon_group,Trainer);
			obj.Meanscore_CellRidgeparams=Meanscore;
			obj.Stdscore_CellRidgeparams=Stdscore;

			if options.SaveWeights
				obj.Weights_ridge=B_group;
			end

			if options.SaveDemoCellActivity
				temp=Meanscore.Rsquared;
				a=sum(temp,1);
				[~,ridgeparam2save]=max(a);
				y=y_group{1};
				y_recon=y_recon_group{1}(:,:,ridgeparam2save);
			end

			if options.saveCellActivity
				obj.Resp_group=y_group;
				obj.Rec_group=y_recon_group;
			end

			if ~isempty(options.Directory2saveData)
				currdir=append(options.Directory2saveData,"_",nameOfCurrentSetting,".mat");
				csp=copy(obj);
				csp.clearData(Properties2keep=["Meanscore_CellRidgeparams","Stdscore_CellRidgeparams","ridgeparams"])
				if options.DisplayTimeTook2Save
					tic
					if ~options.SaveDemoCellActivity
						% save(currdir,"Meanscore","Stdscore");
						save(currdir,"csp");
					else
						save(currdir,"csp","y","y_recon");
					end
					fprintf("Time taken to save data: ")
					toc
				else
					if ~options.SaveDemoCellActivity
						save(currdir,"csp");
					else
						save(currdir,"csp","y","y_recon");
					end
				end
			end
			
			function [currResult,B,y,y_recon] = processGroup(Trainer, Tester, ridgeparams, MakeDataSingle, groupIndex)
				% Process trainer data
				if MakeDataSingle
					y = gpuArray(single(Trainer.Resp{groupIndex}));
					x = gpuArray(single(Trainer.Pred{groupIndex}));
				else
					y = gpuArray(Trainer.Resp{groupIndex});
					x = gpuArray(Trainer.Pred{groupIndex});
				end
				B = ridge_optimized(y, x, ridgeparams);
				clear y x;

				% Process tester data
				if MakeDataSingle
					x = gpuArray(single(Tester.Pred{groupIndex}));
					y = gpuArray(single(Tester.Resp{groupIndex}));
				else
					x = gpuArray(Tester.Pred{groupIndex});
					y = gpuArray(Tester.Resp{groupIndex});
				end
				y_recon = B(1,:,:) + pagemtimes(x, B(2:end,:,:));
				currResult = gather(fastPearsonCorr(y, y_recon, 1));

				B=gather(B);
				y=gather(y);
				y_recon=gather(y_recon);
			end

			function [Meanscore,Stdscore]=calculateRegressionPerformance(y_group,y_recon_group,Trainer)
				y_group=permute(y_group(:),[2,3,4,1]);
				y_recon_group=permute(y_recon_group(:),[2,3,4,1]);
				y_train_group=permute(Trainer.Resp(:),[2,3,4,1]);
				Meanscore=[];
				Stdscore=[];

				[Meanscore.MeanSquaredError,Stdscore.MeanSquaredError]=AssessPerformanceAllCell(y_group,y_recon_group,@MeanSquaredError);
				[Meanscore.MeanAbsoluteError,Stdscore.MeanAbsoluteError]=AssessPerformanceAllCell(y_group,y_recon_group,@MeanAbsoluteError);
				% [Meanscore.Rsquared,Stdscore.Rsquared]=AssessPerformanceAllCell(y_group,y_recon_group,@Rsquared);
				[Meanscore.Rsquared,Stdscore.Rsquared]=AssessRsquared(y_group,y_recon_group,y_train_group);
				[Meanscore.ExplainedVariance,Stdscore.ExplainedVariance]=AssessPerformanceAllCell(y_group,y_recon_group,@ExplainedVariance);
				[Meanscore.PearsonCorrelation,Stdscore.PearsonCorrelation]=AssessPerformanceAllCell(y_group,y_recon_group,@PearsonCorrelation);
			end
			function [meanscore,stdscore]=AssessRsquared(y_group,y_recon_group,y_train_group)
				score_group=cellfun(@Rsquared,y_group,y_recon_group,y_train_group,UniformOutput=false);
				score_group=cell2mat(score_group);	% 1 x cell x ridgeparams x groups
				score_group=DeleteSingleDim(score_group,1);	% cell x ridgeparams x groups
				meanscore=mean(score_group,3);
				stdscore=std(score_group,1,3);
			end
			function [meanscore,stdscore]=AssessPerformanceAllCell(y_group,y_recon_group,func_handle)
				score_group=cellfun(func_handle,y_group,y_recon_group,UniformOutput=false);
				score_group=cell2mat(score_group);	% 1 x cell x ridgeparams x groups
				score_group=DeleteSingleDim(score_group,1);	% cell x ridgeparams x groups
				meanscore=mean(score_group,3);
				stdscore=std(score_group,1,3);
			end
			% Below, first dim is time, second is cell, third is
			% ridgeparams
			% --> result: 1 x cell x ridgeparams
			function score = MeanSquaredError(y, y_recon)
				score = sum((y - y_recon).^2, 1) / size(y, 1);
			end
			function score = MeanAbsoluteError(y, y_recon)
				score = sum(abs(y - y_recon), 1) / size(y, 1);
			end
			%%%%%%%%%%%
			% function score = Rsquared(y, y_recon)
			% 	y = y + zeros(size(y_recon), class(y_recon));  % broadcast to match dimensions
			% 	normy = y - mean(y, 1);
			% 	score = 1 - ( sum((y - y_recon).^2, 1) ./ sum(normy.^2, 1) );
			% end
			%%%%%%%%%%%
			function score = Rsquared(y, y_recon,y_trainer)
				y = y + zeros(size(y_recon), class(y_recon));  % broadcast to match dimensions
				y_trainer=cast(y_trainer,class(y_recon));
				mean_y=mean(y_trainer,1);
				normy = y - mean_y;
				score = 1 - ( sum((y - y_recon).^2, 1) ./ sum(normy.^2, 1) );
			end
			function score = ExplainedVariance(y, y_recon)
				y = y + zeros(size(y_recon), class(y_recon));  % broadcasting
				score=1-(var(y-y_recon,0,1)./var(y,0,1));
				% score = 1 - ( sum((y - y_recon).^2, 1) ./ sum(y.^2, 1) );
			end
			function score = PearsonCorrelation(y, y_recon)
				score = fastPearsonCorr(y, y_recon,1);
			end
		end
	end
	properties
		Weights_ridge_mean
	end
	methods
		function take_mean_of_weights_Ridge(obj)
			w=obj.Weights_ridge;
			w=permute(w(:),[2,3,4,1]);
			w=cell2mat(w);
			meanw=mean(w,4);
			obj.Weights_ridge_mean=meanw;
		end
	end
	%%%%%%%
	%% new ridge regression, compatible with ablation and weight zeroing
	properties
		BodyPart_4ridge
		Ridge_y
		Ridge_y_rec
	end
	methods
		function ridgeResults=Ridge_Train(obj,ridgeparams,nameOfCurrentSetting,options)
			arguments
				obj
				ridgeparams
				nameOfCurrentSetting (1,1) string
				options.MakeDataSingle=false
				options.DisplayGPUMemory=false
				options.SaveWeights	(1,1) logical =false
				options.indices_of_predictors_to_use	=[]
				options.indices_to_zeroOutWeights=[]
				options.useGPU (1,1) =true
			end
			obj.nameOfCurrentSettingRidge=nameOfCurrentSetting;
			Trainer=obj.Trainer;
			Tester=obj.Tester;
			B_group=[];
			y_group=[];
			y_recon_group=[];
			y_null_group=[];
			%%%%%%%%%%%%%%%%%%
			for i=1:obj.groupnum
				if options.DisplayGPUMemory
					fprintf("Available GPU Memory: %d\n", gpuDevice().AvailableMemory);
				end
				 [B_group{i},y_group{i},y_recon_group{i},y_null_group{i}] = processGroup(Trainer, Tester, ridgeparams, options.MakeDataSingle, i,options);
			end
			obj.ridgeparams=ridgeparams;

			Meanscore=assessPerformance(y_group,y_recon_group,y_null_group);
			obj.Meanscore_CellRidgeparams=Meanscore;
			%%%%%%%%%%%%%%%%%%
			if options.SaveWeights
				obj.Weights_ridge=B_group;
			end

			obj.Ridge_y=y_group;
			obj.Ridge_y_rec=y_recon_group;
			
			function [B,y,y_recon,y_null] = processGroup(Trainer, Tester, ridgeparams, MakeDataSingle, groupIndex,options)
				% Process trainer data
				if MakeDataSingle
					y = (single(Trainer.Resp{groupIndex}));
					x = (single(Trainer.Pred{groupIndex}));
				else
					y = (Trainer.Resp{groupIndex});
					x = (Trainer.Pred{groupIndex});
				end
				if options.useGPU
					y=gpuArray(y);
					x=gpuArray(x);
				end
				if ~isempty(options.indices_of_predictors_to_use)
					xind=makeIND(x,2,options.indices_of_predictors_to_use);
					x=x(xind{:});
					if obj.BodyPart_4ridge=="All"
						error("when indices_of_predictors_to_use is specified, obj.BodyPart_4ridge is not supposed to be the name of the body part, not All")
					end
				else
					obj.BodyPart_4ridge="All";
				end
				B = ridge_optimized(y, x, ridgeparams);
				y_train=gather(y);
				clear y x;

				% Process tester data
				if MakeDataSingle
					x = (single(Tester.Pred{groupIndex}));
					y = (single(Tester.Resp{groupIndex}));
				else
					x = (Tester.Pred{groupIndex});
					y = (Tester.Resp{groupIndex});
				end
				if options.useGPU
					y=gpuArray(y);
					x=gpuArray(x);
				end
				if ~isempty(options.indices_of_predictors_to_use)
					xind=makeIND(x,2,options.indices_of_predictors_to_use);
					x=x(xind{:});
				end

				y_recon = B(1,:,:) + pagemtimes(x, B(2:end,:,:));

				B=gather(B);
				y=gather(y);
				y_recon=gather(y_recon);				%% testing time points x cell x ridgeparams
				y=y+zeros(size(y_recon));				%% testing time points x cell x ridgeparams
				y_null=make_null_model(y_recon,y_train);%% testing time points x cell x ridgeparams
			end

			function y_train_mean=make_null_model(y_recon,y_train)
				y_train_mean=mean(y_train,1);					%% training time points x cell
				y_train_mean=y_train_mean+zeros(size(y_recon));	%% testing time points x cell x ridgeparams
			end

			function Meanscore=assessPerformance(y_group,y_recon_group,y_null_group)
				% arrays inside are cell x ridgeparams
				y=cell2mat(y_group(:));
				y_recon=cell2mat(y_recon_group(:));
				y_null=cell2mat(y_null_group(:));
				Meanscore.MeanSquaredError=MeanSquaredError(y, y_recon);
				Meanscore.MeanAbsoluteError=MeanAbsoluteError(y, y_recon);
				Meanscore.Rsquared=Rsquared(y, y_recon);
				Meanscore.pseudoRsquared=pseudoRsquared(y, y_recon,y_null);
				Meanscore.ExplainedVariance=ExplainedVariance(y, y_recon);
				Meanscore.PearsonCorrelation=PearsonCorrelation(y, y_recon);
			end

			
			% Below, first dim is time, second is cell, third is
			% ridgeparams
			% --> result: cell x ridgeparams
			function score = MeanSquaredError(y, y_recon)
				score = sum((y - y_recon).^2, 1) / size(y, 1);
				score=DeleteSingleDim(score,1);
			end
			function score = MeanAbsoluteError(y, y_recon)
				score = sum(abs(y - y_recon), 1) / size(y, 1);
				score=DeleteSingleDim(score,1);
			end
			function score = Rsquared(y, y_recon)
				normy = y - mean(y, 1);
				score = 1 - ( sum((y - y_recon).^2, 1) ./ sum(normy.^2, 1) );
				score=DeleteSingleDim(score,1);
			end
			function score = pseudoRsquared(y, y_recon,y_null)
				score=1-(sum((y-y_recon).^2,1)./sum((y-y_null).^2,1));	% 1 x cell x ridge params
				score=DeleteSingleDim(score,1);
			end
			function score = ExplainedVariance(y, y_recon)
				score=1-(var(y-y_recon,0,1)./var(y,0,1));
				score=DeleteSingleDim(score,1);
			end
			function score = PearsonCorrelation(y, y_recon)
				score = fastPearsonCorr(y, y_recon,1);
				score=DeleteSingleDim(score,1);
			end
		end
	end
	%% make_body_part_indices_for_ridge %%%%%%
	properties
		Body_part_indices_for_ridge
	end
	methods
		function all_logs=make_body_part_indices_for_ridge(obj,options)
			arguments
				obj
				options.names=[]
				options.BodyPart (1,1) {mustBeMember(options.BodyPart,[ ...
											"RightHand", "RightElbow", ...
											"RightToe", "RightHeel", "RightKnee", ...
											"BackU", "BackM", "BackL", "Tail", ...
											"LeftHand", "LeftElbow", "LeftToe", ...
											"LeftHeel", "LeftKnee","All"])}	="All"
				options.XYZ_ANG_DOT (1,1){mustBeMember(options.XYZ_ANG_DOT,["XYZ","ANG","DOT","All"])}="All"
				options.XYZ_coord (1,1) {mustBeMember(options.XYZ_coord,["x","y","z","All"])}="All"
				options.Derivatives (1,1) string {mustBeMember(options.Derivatives,["0","1","2","All"])}="All"
				% options.speed_and_accel (1,1) {mustBeMember(options.speed_and_accel,["speed","accel","speed and accel","None"])}="speed and accel"
			end
			obj.BodyPart_4ridge=options.BodyPart;
			if isempty(options.names)
				all_names=obj.all_GaussianHot_names(:)';
			else
				all_names=options.names;
			end
			all_logs={};
			all_logs{end+1}=asessBodyPart(all_names,options);
			all_logs{end+1}=assessXYZ_ANG_DOT(all_names,options);
			all_logs{end+1}=assessXYZ_coord(all_names,options);
			all_logs{end+1}=assessDerivatives(all_names,options);
			all_logs=cell2mat(all_logs(:));
			all_logs=all(all_logs,1);
			obj.Body_part_indices_for_ridge=all_logs;
			if ~any(all_logs,"all")
				warning("Returning empty indices...")
			end
			

			function currlog=asessBodyPart(all_names,options)
				if options.BodyPart=="All"
					currlog=true(size(all_names));
					return;
				end
				currlog=containsSubstring(all_names,options.BodyPart);				
			end
			function currlog=assessXYZ_ANG_DOT(all_names,options)
				if options.XYZ_ANG_DOT=="All"
					currlog=true(size(all_names));
					return;
				end
				currlog=containsSubstring(all_names,options.XYZ_ANG_DOT);
			end
			function currlog=assessXYZ_coord(all_names,options)
				if options.XYZ_coord=="All"
					currlog=true(size(all_names));
					return;
				end
				currlog=containsSubstring(all_names,"_"+options.XYZ_coord,CapSensitive=true);
			end
			function currlog=assessDerivatives(all_names,options)
				if options.Derivatives=="All"
					currlog=true(size(all_names));
					return;
				end
				positstr=["tXYZbp","tDOTbp","tANGbp"];
				speedstr=["tDXYZbp","tDDOTbp","tDANGbp"];
				accelstr=["tD2XYZbp","tD2DOTbp","tD2ANGbp"];
				switch options.Derivatives
					case "0"
						currlog=containsSubstring(all_names,positstr);
					case "1"
						currlog=containsSubstring(all_names,speedstr);
					case "2"
						currlog=containsSubstring(all_names,accelstr);
				end
			end
		end

		function all_logs=search_body_part_indices(obj,name_list,options)
			arguments
				obj
				name_list
				options.BodyPart (1,1) {mustBeMember(options.BodyPart,[ ...
											"RightHand", "RightElbow", ...
											"RightToe", "RightHeel", "RightKnee", ...
											"BackU", "BackM", "BackL", "Tail", ...
											"LeftHand", "LeftElbow", "LeftToe", ...
											"LeftHeel", "LeftKnee","All"])}	="All"
				options.XYZ_ANG_DOT (1,1){mustBeMember(options.XYZ_ANG_DOT,["XYZ","ANG","DOT","All"])}="All"
				options.XYZ_coord (1,1) {mustBeMember(options.XYZ_coord,["x","y","z","All"])}="All"
				options.Derivatives (1,1) string {mustBeMember(options.Derivatives,["0","1","2","All"])}="All"
				% options.speed_and_accel (1,1) {mustBeMember(options.speed_and_accel,["speed","accel","speed and accel","None"])}="speed and accel"
			end
			all_names=name_list(:);
			all_logs={};
			all_logs{end+1}=asessBodyPart(all_names,options);
			all_logs{end+1}=assessXYZ_ANG_DOT(all_names,options);
			all_logs{end+1}=assessXYZ_coord(all_names,options);
			all_logs{end+1}=assessDerivatives(all_names,options);
			all_logs=cat(2,all_logs{:});
			all_logs=all(all_logs,2);
			all_logs=reshape(all_logs,size(name_list));
			if ~any(all_logs,"all")
				warning("Returning empty indices...")
			end
			function currlog=asessBodyPart(all_names,options)
				if options.BodyPart=="All"
					currlog=true(size(all_names));
					return;
				end
				currlog=containsSubstring(all_names,options.BodyPart);				
			end
			function currlog=assessXYZ_ANG_DOT(all_names,options)
				if options.XYZ_ANG_DOT=="All"
					currlog=true(size(all_names));
					return;
				end
				currlog=containsSubstring(all_names,options.XYZ_ANG_DOT);
			end
			function currlog=assessXYZ_coord(all_names,options)
				if options.XYZ_coord=="All"
					currlog=true(size(all_names));
					return;
				end
				currlog=containsSubstring(all_names,"_"+options.XYZ_coord,CapSensitive=true);
			end
			function currlog=assessDerivatives(all_names,options)
				if options.Derivatives=="All"
					currlog=true(size(all_names));
					return;
				end
				positstr=["tXYZbp","tDOTbp","tANGbp"];
				speedstr=["tDXYZbp","tDDOTbp","tDANGbp"];
				accelstr=["tD2XYZbp","tD2DOTbp","tD2ANGbp"];
				switch options.Derivatives
					case "0"
						currlog=containsSubstring(all_names,positstr);
					case "1"
						currlog=containsSubstring(all_names,speedstr);
					case "2"
						currlog=containsSubstring(all_names,accelstr);
				end
			end
		end
	end






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	%%%% Analysis of ridge performance
% 	methods
% 		function compareLength_warningUponMismatch(obj,input1,input2)
% 			if length(input1)~=length(input2)
% 				warning("Input length does not match")
% 			end
% 		end
% 		function plotOptimalRidgeParam(obj,ridgeSetting2plot)
% 			arguments
% 				obj
% 				ridgeSetting2plot (1,1) string
% 			end
% 			ridgeResults_cell=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot);
% 			ridgeparams=obj.ridgeparams;
% 
% 			if length(ridgeResults_cell)~=length(obj.info.NoAltplane.nameFull_p)
% 				warning("ridgeResults_cell do not contain the same number of elements as the number of planes")
% 			end
% 
% 			ax=gobjects(0);
% 			for i=1:length(ridgeResults_cell)
% 				ridgeResults=ridgeResults_cell{i};
% 				% For each cell, find the index of the ridge parameter that maximizes the correlation
% 				[~, maxIdx] = max(ridgeResults, [], 2);
% 
% 				% Count the frequency for each parameter index (assuming indices 1:numel(ridgeparams))
% 				counts = accumarray(maxIdx, 1, [numel(ridgeparams) 1]);
% 
% 				% Plot the histogram as a line plot with x-axis in log scale
% 				ax(i)=nexttile;
% 				semilogx(ridgeparams, counts, '-o', 'LineWidth', 2);
% 				xlabel('Ridge Parameter Value (log scale)');
% 				ylabel('Number of Cells');
% 				grid on;
% 				title(obj.info.NoAltplane.nameFull_p(i));
% 			end
% 			linkaxes(ax,"x");
% 			Fig=gcf;
% 			Fig.Children.Title.String=append('Frequency Distribution of Ridge Regression Parameter Maximizing Correlation for model: ',regexprep(ridgeSetting2plot,'_',' '));
% 		end
% 		%%
% 		function plotOptimalRidgeParam_ShiftTimePointsOfCellData(obj,ridgeSetting2plot)
% 			arguments
% 				obj
% 				ridgeSetting2plot
% 			end
% 			ridgedata=obj.ridgeResults_setting_plane_CellRidgeShift.(ridgeSetting2plot);
% 			counter=[];
% 			for i=1:length(ridgedata)
% 				data=ridgedata{i};
% 				[~,maxind]=max(data,[],2);
% 				maxIdx=squeeze(maxind);
% 
% 				% Get dimensions
% 				[numCells, numConds] = size(maxIdx);
% 				% Create condition indices for each cell (each column corresponds to a condition)
% 				condIndices = repmat(1:numConds, numCells, 1);
% 				% Flatten the matrices so that each row corresponds to one observation
% 				maxIdx_flat = maxIdx(:);         % Ridge parameter index per cell-condition
% 				condIndices_flat = condIndices(:); % Corresponding condition index
% 				% Use accumarray to compute counts:
% 				% Rows correspond to ridge parameter indices, columns to conditions.
% 				counter{i} = accumarray([maxIdx_flat, condIndices_flat], 1, [numel(obj.ridgeparams), numConds]);
% 			end
% 			Color=parula(numel(obj.shiftIn10Hz_vec));
% 			Color=num2cell(Color,2);
% 
% 			obj.plotTile_allPlane(@plotfunc,counter,obj.info.NoAltplane.nameFull_p, ...
% 				Legend=arrayfun(@(c) sprintf('Time Shift %g sec', c), obj.shiftIn10Hz_vec/10, 'UniformOutput', false),LegendLoc="east");
% 
% 			function plotfunc(counter,titletext)
% 				h=semilogx(obj.ridgeparams, counter, '-', 'LineWidth', 2);
% 				set(h,{'Color'},Color)
% 				xlabel('Ridge Parameter Value (log scale)');
% 				ylabel('Number of Cells');
% 				title('Distribution of Optimal Ridge Parameters Across Conditions');
% 				grid on;
% 				title(titletext);
% 			end
% 		end
% 		%%
% 		function plotAvgFittingPerformance(obj,ridgeSetting2plot,options)
% 			arguments
% 				obj
% 				ridgeSetting2plot (1,1) string
% 				options.linkaxes (1,1) string {mustBeMember(options.linkaxes,["x","y","xy"])}=[]
% 				options.FontSize=[]
% 			end
% 			options=namedargs2cell(options);
% 			ridgeResults_cell=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot);
% 			ridgeparams=obj.ridgeparams;
% 			obj.compareLength_warningUponMismatch(ridgeResults_cell,obj.info.NoAltplane.nameFull_p);
% 			obj.plotTile_allPlane(@plotfunc,ridgeResults_cell,obj.info.NoAltplane.nameFull_p,options{:})
% 			Fig=gcf;
% 			Fig.Children.Title.String=append('Average Fitting Performance for Each Ridge Parameter for model: ',regexprep(ridgeSetting2plot,'_',' '));
% 			function plotfunc(Data,titletext)
% 				% Compute the mean correlation across cells for each parameter.
% 				avgCorr = mean(Data, 1);
% 				% Create a figure and plot using a logarithmic x-axis.
% 				semilogx(ridgeparams, avgCorr, '-o', 'LineWidth', 2);
% 				% Label the axes and add a title.
% 				xlabel('Ridge Parameter (log scale)');
% 				ylabel('Average Pearson Correlation');
% 				grid on;
% 				title(titletext);
% 			end
% 		end
% 		%%
% 		function plotAvgFittingPerformance_ShiftTimePointsOfCellData(obj,ridgeSetting2plot,options)
% 			arguments
% 				obj
% 				ridgeSetting2plot (1,1) string
% 				options.linkaxes (1,1) string {mustBeMember(options.linkaxes,["x","y","xy"])}=[]
% 				options.FontSize=[]
% 			end
% 			options=namedargs2cell(options);
% 			ridgeResuls=obj.ridgeResults_setting_plane_CellRidgeShift.(ridgeSetting2plot);
% 			Color=parula(numel(obj.shiftIn10Hz_vec));
% 			Color=num2cell(Color,2);
% 			obj.compareLength_warningUponMismatch(ridgeResuls,obj.info.NoAltplane.nameFull_p);
% 			obj.plotTile_allPlane(@plotfunc,ridgeResuls,obj.info.NoAltplane.nameFull_p,options{:}, ...
% 				Legend=arrayfun(@(c) sprintf('Time Shift %g sec', c), obj.shiftIn10Hz_vec/10, 'UniformOutput', false),LegendLoc="east");
% 			Fig=gcf;
% 			Fig.Children.Title.String=append('Average Fitting Performance for Each Ridge Parameter for model: ',regexprep(ridgeSetting2plot,'_',' '));
% 
% 			function plotfunc(Data,titletext)
% 				% Compute the mean correlation across cells for each parameter.
% 				avgCorr = mean(Data, 1);	% 1 x ridgeparams x shift
% 				avgCorr=squeeze(avgCorr);
% 				% Create a figure and plot using a logarithmic x-axis.
% 				h=semilogx(obj.ridgeparams, avgCorr, '-', 'LineWidth', 2);
% 				set(h,{'Color'},Color);
% 				% Label the axes and add a title.
% 				xlabel('Ridge Parameter (log scale)');
% 				ylabel('Average Pearson Correlation');
% 				grid on;
% 				title(titletext);
% 			end
% 		end
% 		%%
% 		function plotCorrelationDistribution(obj,ridgeSetting2plot,ridgeparam2plot,options)
% 			arguments
% 				obj
% 				ridgeSetting2plot (1,1) string
% 				ridgeparam2plot=[]
% 				options.linkaxes (1,1) string {mustBeMember(options.linkaxes,["x","y","xy"])}=[]
% 				options.FontSize=[]
% 			end
% 			options=namedargs2cell(options);
% 			ridgeResults_cell=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot);
% 			ridgeparams=obj.ridgeparams;
% 			obj.compareLength_warningUponMismatch(ridgeResults_cell,obj.info.NoAltplane.nameFull_p);
% 			obj.plotTile_allPlane(@plotfunc,ridgeResults_cell,obj.info.NoAltplane.nameFull_p,options{:})
% 			Fig=gcf;
% 			Fig.Children.Title.String=append("Distribution of correlation for model:  ",regexprep(ridgeSetting2plot,'_',' '));
% 
% 			function plotfunc(ridgeResults,titletext)
% 				% Compute the average correlation across cells for each ridge parameter.
% 				avgCorr = mean(ridgeResults, 1);
% 
% 				% Identify the ridge parameter that maximizes the average correlation.
% 				if isempty(ridgeparam2plot)
% 					[maxAvgCorr, bestIdx] = max(avgCorr);
% 					bestParam = ridgeparams(bestIdx);
% 				else
% 					bestParam=ridgeparam2plot;
% 					bestIdx=find(ridgeparams==bestParam,1);
% 				end
% 
% 				% Subplot 2: Distribution of correlation coefficients for the best parameter.
% 				% Get the distribution (across cells) for the best ridge parameter.
% 				bestDistribution = ridgeResults(:, bestIdx);
% 				histogram(bestDistribution, 'Normalization', 'pdf');
% 				xlabel('Pearson Correlation');
% 				ylabel('Probability Density');
% 				if isempty(ridgeparam2plot)
% 					% title(sprintf('Distribution for Best Parameter: %.2g', bestParam));
% 					text(0.4, 0.95, sprintf('Best Parameter: %.2g', bestParam), ...
% 						'Units', 'normalized', 'VerticalAlignment', 'top', ...
% 						'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight', 'bold', ...
% 						'BackgroundColor', 'white', 'EdgeColor', 'black');
% 				else
% 					% title(sprintf('Distribution current Parameter: %.2g', bestParam));
% 					text(0.4, 0.95, sprintf('Current Parameter: %.2g', bestParam), ...
% 						'Units', 'normalized', 'VerticalAlignment', 'top', ...
% 						'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight', 'bold', ...
% 						'BackgroundColor', 'white', 'EdgeColor', 'black');
% 				end
% 				grid on;hold on;
% 				xline(0,'Color','red',LineWidth=2);hold off;
% 				title(sprintf('ridge param: %.2g, %s',bestParam,titletext));
% 			end
% 		end
% 		%%
% 		function plotCorrelationDistribution_ShiftTimePointsOfCellData(obj,ridgeSetting2plot,ridgeparam2plot,options)
% 			arguments
% 				obj
% 				ridgeSetting2plot (1,1) string
% 				ridgeparam2plot=[]
% 				options.linkaxes (1,1) string {mustBeMember(options.linkaxes,["","x","y","xy"])}=""
% 				options.FontSize=[]
% 			end
% 			options=namedargs2cell(options);
% 			ridgeResults_cell=obj.ridgeResults_setting_plane_CellRidgeShift.(ridgeSetting2plot);
% 			ridgeparams=obj.ridgeparams;
% 			Color=parula(numel(obj.shiftIn10Hz_vec));
% 			Color=num2cell(Color,2);
% 			obj.compareLength_warningUponMismatch(ridgeResults_cell,obj.info.NoAltplane.nameFull_p);
% 			obj.plotTile_allPlane(@plotfunc,ridgeResults_cell,obj.info.NoAltplane.nameFull_p,options{:})
% 			ChangeFigureSettings(MasterTitle=append("Distribution of correlation for model:  ",regexprep(ridgeSetting2plot,'_',' ')), ...
% 				Legend=arrayfun(@(c) sprintf('Time Shift %g sec', c), obj.shiftIn10Hz_vec/10, 'UniformOutput', false),LegendLoc="east")
% 			% Fig.Children.Title.String=append("Distribution of correlation for model:  ",regexprep(ridgeSetting2plot,'_',' '));
% 
% 			function plotfunc(ridgeResults,titletext)
% 				data=[];	% cell x shift
% 				for i=1:size(ridgeResults,3)
% 					currdata=ridgeResults(:,:,i);
% 
% 					% Identify the ridge parameter that maximizes the average correlation.
% 					if isempty(ridgeparam2plot)
% 						% Compute the average correlation across cells for each ridge parameter.
% 						avgCorr = mean(currdata, 1);
% 						[maxAvgCorr, bestIdx] = max(avgCorr);
% 						bestParam = ridgeparams(bestIdx);
% 					else
% 						bestParam=ridgeparam2plot;
% 						bestIdx=find(ridgeparams==bestParam,1);
% 					end
% 
% 					% Subplot 2: Distribution of correlation coefficients for the best parameter.
% 					% Get the distribution (across cells) for the best ridge parameter.
% 					bestDistribution = currdata(:, bestIdx);
% 					data(:,i)=bestDistribution;
% 				end
% 				[~,edges]=histcounts(data);
% 				pdata=[];
% 				for i=1:size(ridgeResults,3)
% 					pdata(:,i)=histcounts(data(:,i),edges);
% 				end
% 				x=edges(1:end-1)+diff(edges,1);
% 				h=plot(x(:), pdata, '-', 'LineWidth', 1);
% 				set(h,{'Color'},Color);
% 				xlabel('Pearson Correlation');
% 				ylabel('Probability Density');
% 				if isempty(ridgeparam2plot)
% 					% title(sprintf('Distribution for Best Parameter: %.2g', bestParam));
% 					text(0.4, 0.95, sprintf('Best Parameter: %.2g', bestParam), ...
% 						'Units', 'normalized', 'VerticalAlignment', 'top', ...
% 						'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight', 'bold', ...
% 						'BackgroundColor', 'white', 'EdgeColor', 'black');
% 				else
% 					% title(sprintf('Distribution current Parameter: %.2g', bestParam));
% 					text(0.4, 0.95, sprintf('Current Parameter: %.2g', bestParam), ...
% 						'Units', 'normalized', 'VerticalAlignment', 'top', ...
% 						'HorizontalAlignment', 'left', 'FontSize', 8, 'FontWeight', 'bold', ...
% 						'BackgroundColor', 'white', 'EdgeColor', 'black');
% 				end
% 				grid on;hold on;
% 				xline(0,'Color','red',LineWidth=2);hold off;
% 				title(sprintf('ridge param: %.2g, %s',bestParam,titletext));
% 			end
% 		end
% %%
% 		function plotCompareRidge(obj,ridgeSetting2plot1,ridgeSetting2plot2,ridgeparam2plot,plotoptions,options)
% 			arguments
% 				obj
% 				ridgeSetting2plot1 (1,1) string
% 				ridgeSetting2plot2 (1,1) string
% 				ridgeparam2plot=[]
% 				plotoptions.SkewnessThresh=[]
% 				options.linkaxes (1,1) string {mustBeMember(options.linkaxes,["","x","y","xy"])}=""
% 				options.FontSize=[]
% 			end
% 			SkewnessThresh=plotoptions.SkewnessThresh;
% 			options=namedargs2cell(options);
% 			ridgeResults1=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot1);
% 			ridgeResults2=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot2);
% 			ridgeparams=obj.ridgeparams;
% 			text1=regexprep(ridgeSetting2plot1,'_',' ');
% 			text2=regexprep(ridgeSetting2plot2,'_',' ');
% 			Skewlog=generateSkewList();
% 
% 			obj.plotTile_allPlane(@plotfunc,ridgeResults1,ridgeResults2,obj.info.NoAltplane.nameFull_p,Skewlog,options{:});
% 			axis equal
% 			Fig=gcf;
% 			Fig.Children.Title.String=sprintf("Correlation coefficient comparison for models:  [%s] VS [%s], Skewness>%g",regexprep(ridgeSetting2plot1,'_',' '),regexprep(ridgeSetting2plot2,'_',' '),SkewnessThresh);
% 
% 			function Skewlog=generateSkewList()
% 				Skewlog=[];
% 				for i=1:length(ridgeResults1)
% 					obj.plane=i;
% 					obj.load_specified(obj.plane,"Skew");
% 					% filter out cells with low skewness
% 					if ~isempty(SkewnessThresh)
% 						Skewlog{i}=obj.Skew(:)>=SkewnessThresh;
% 					else
% 						Skewlog{i}=logical(ones(size(obj.Skew(:))));
% 					end
% 				end
% 			end
% 
% 			function plotfunc(ridgeResults1,ridgeResults2,titletext,skewlog)
% 				IND1=extractBestRidgeParam(ridgeResults1);
% 				IND2=extractBestRidgeParam(ridgeResults2);
% 				if ~isempty(ridgeparam2plot)
% 					data1=ridgeResults1(skewlog,IND1);
% 					data2=ridgeResults2(skewlog,IND2);
% 					scatter(ridgeResults1(:,IND1),ridgeResults2(:,IND2),5,"filled");
% 					xlabel(sprintf('%s, ridge parameter: %.2g', text1,ridgeparams(IND1)))
% 					ylabel(sprintf('%s, ridge parameter: %.2g', text2,ridgeparams(IND2)))
% 				else
% 					data1=ridgeResults1(sub2ind(size(ridgeResults1),[1:size(ridgeResults1)]',IND1(:)));
% 					data2=ridgeResults2(sub2ind(size(ridgeResults2),[1:size(ridgeResults2)]',IND2(:)));
% 					data1=data1(skewlog,:);
% 					data2=data2(skewlog,:);
% 					scatter(ridgeResults1(sub2ind(size(ridgeResults1),[1:size(ridgeResults1)]',IND1(:))),...
% 						ridgeResults2(sub2ind(size(ridgeResults2),[1:size(ridgeResults2)]',IND2(:))),5,"filled");
% 				end
% 				xlabel(append(text1))
% 				ylabel(append(text2))
% 				hold on;plot([-0.2,1],[-0.2,1],Color='red',LineWidth=1);
% 				% Fit a linear model
% 				mdl = fitlm(data1, data2);
% 				% Extract coefficients: Intercept and slope
% 				coeffs = mdl.Coefficients.Estimate;
% 				intercept = coeffs(1);
% 				slope = coeffs(2);
% 				% Extract the p-value for the slope coefficient
% 				p_val = mdl.Coefficients.pValue(2);
% 				% Generate fitted line data
% 				x_fit = linspace(-0.2, 1, 100)';
% 				y_fit = intercept + slope*x_fit;
% 
% 				% Plot the fitted line
% 				plot(x_fit, y_fit, 'blue', 'LineWidth', 2);
% 
% 				% Add legend with the p-value (formatted to 3 decimal places)
% 				eqStr = sprintf('Fit: y = %.2fx + %.2f, p = %.3f', slope, intercept, p_val);
% 				legend('Data','Partition', eqStr, 'Location', 'best');
% 
% 				hold off;
% 				title(sprintf("%s, cell count: %d / %d",titletext,sum(skewlog),length(skewlog)))
% 
% 				function ridgeIND=extractBestRidgeParam(ridgeResults)
% 					% Identify the ridge parameter that maximizes the average correlation.
% 					if isempty(ridgeparam2plot)
% 						[maxAvgCorr, bestIdx] = max(ridgeResults,[],2);
% 						bestParam = ridgeparams(bestIdx);
% 					else
% 						bestParam=ridgeparam2plot;
% 						bestIdx=find(ridgeparams==bestParam,1);
% 					end
% 					ridgeIND=bestIdx;
% 				end
% 			end
% 		end
% 		%%
% 		function plotCompareRidge_ShiftTimePointsOfCellData(obj,ridgeSetting2plot1,ridgeSetting2plot2,ridgeparam2plot,shiftVal2plotInSec,plotoptions,options)
% 			arguments
% 				obj
% 				ridgeSetting2plot1 (1,1) string
% 				ridgeSetting2plot2 (1,1) string
% 				ridgeparam2plot=[]
% 				shiftVal2plotInSec=[]
% 
% 				plotoptions.SkewnessThresh=[]
% 				options.linkaxes {mustBeMember(options.linkaxes,["","x","y","xy"])}=""
% 				options.FontSize=[]
% 			end
% 			SkewnessThresh=plotoptions.SkewnessThresh;
% 			shiftind=find(obj.shiftIn10Hz_vec==shiftVal2plotInSec*10);
% 
% 			csp=copy(obj);
% 			csp.ridgeResults_setting_plane_CellRidge=applyFunctionLayers(obj.ridgeResults_setting_plane_CellRidgeShift,@(a)a(:,:,shiftind));
% 			csp.plotCompareRidge(ridgeSetting2plot1,ridgeSetting2plot2,ridgeparam2plot,SkewnessThresh=SkewnessThresh);
% 			ChangeFigureSettings(MasterTitle=sprintf("Correlation coefficient comparison for models:  [%s] VS [%s], Skewness>%g, Time shift: %g sec", ...
% 				regexprep(ridgeSetting2plot1,'_',' '),regexprep(ridgeSetting2plot2,'_',' '),SkewnessThresh,shiftVal2plotInSec), ...
% 				Linkaxes="xyz",MasterTitleSize=20);
% 		end
% 		%%
% 		function plotRidgeChangeAcrossShift_ShiftTimePointsOfCellData(obj,ridgeSetting2plot,ridgeparam2plot,options)
% 			arguments
% 				obj
% 				ridgeSetting2plot (1,1) string
% 				ridgeparam2plot=[]
% 				options.SkewThresh=[]
% 				options.CorrThresh=[]
% 			end
% 			currdata=obj.ridgeResults_setting_plane_CellRidgeShift.(ridgeSetting2plot);
% 			obj.plotTile_allPlane(@plotfunc,currdata,obj.info.NoAltplane.nameFull_p,obj.generateSkewList(options.SkewThresh));
% 			ChangeFigureSettings(MasterTitle=sprintf("Correlation coefficient over shifting data for model %s, Skew>%g, maxCorr>%g", ...
% 				regexprep(ridgeSetting2plot,"_"," "),options.SkewThresh,options.CorrThresh),Xlabel="Shift (sec)",Ylabel="Correlation coefficient", ...
% 				Linkaxes="x",TitleSize=12);
% 			function plotfunc(currdata,names,Skewlog)
% 				if isempty(ridgeparam2plot)
% 					data=max(currdata,[],2);
% 				else
% 					IND=createIND(currdata,2,find(obj.ridgeparams==ridgeparam2plot));
% 					data=currdata(IND{:});
% 				end
% 				data=DeleteSingleDim(data,2);
% 
% 				sLog=Skewlog;
% 				% if ~isempty(options.SkewThresh)
% 				% 	sLog=Skew>=options.SkewThresh;
% 				% else
% 				% 	sLog=ones(size(data,1),1,"logical");
% 				% end
% 				if ~isempty(options.CorrThresh)
% 					temp=max(data,[],2);
% 					cLog=temp>=options.CorrThresh;
% 				else
% 					cLog=ones(size(data,1),1,"logical");
% 				end
% 				Log=cLog(:).*sLog(:);
% 				data=data(logical(Log),:);
% 
% 
% 				data=data';	% shift x cell
% 				plot(obj.shiftIn10Hz_vec/10,data,Color=[1,1,1]*0.7);
% 				hold on;
% 				plotMeanSE(data,1,x=obj.shiftIn10Hz_vec/10,Color="blue",ShadeAlpha=0.5);
% 				titletext=sprintf("%s, # of cells: %g (%.2g%s)",names,sum(Log),100*sum(Log)/length(Log),"%");
% 				title(titletext);
% 			end
% 		end
% 		%%
% 		function plotRidgeChangeAcrossShiftAll_ShiftTimePointsOfCellData(obj,ridgeparam2plot,options)
% 			arguments
% 				obj
% 				ridgeparam2plot=[]
% 				options.SkewThresh=[]
% 				options.CorrThresh=[]
% 			end
% 			flds=string(fieldnames(obj.ridgeResults_setting_plane_CellRidgeShift));
% 			alldata=struct2cell(obj.ridgeResults_setting_plane_CellRidgeShift);
% 			obj.plotTile_allPlane(@plotfunc,alldata,regexprep(flds,"_"," "));
% 			ChangeFigureSettings(XLabel="Time shift (sec)",Ylabel="Correlation",XLabelSize=20,YlabelSize=20, ...
% 				TitleSize=20,MasterTitle=sprintf("Correlation coefficient over shifting data, Skew>%g, maxCorr>%g, All mice combined", ...
% 				options.SkewThresh,options.CorrThresh),MasterTitleSize=20);
% 
% 			function plotfunc(cdata,titletext)
% 				Skewlog=generateSkewList(obj,options.SkewThresh);
% 				catdata=[];
% 				totalLog=[];
% 				for i=1:length(cdata)
% 					currdata=cdata{i};
% 					if isempty(ridgeparam2plot)
% 						data=max(currdata,[],2);
% 					else
% 						IND=createIND(currdata,2,find(obj.ridgeparams==ridgeparam2plot));
% 						data=currdata(IND{:});
% 					end
% 					data=DeleteSingleDim(data,2);
% 					sLog=Skewlog{i};
% 					if ~isempty(options.CorrThresh)
% 						temp=max(data,[],2);
% 						cLog=temp>=options.CorrThresh;
% 					else
% 						cLog=ones(size(data,1),1,"logical");
% 					end
% 					Log=cLog(:).*sLog(:);
% 					data=data(logical(Log),:);
% 					catdata=cat(1,catdata,data);
% 					totalLog=cat(1,totalLog,Log);
% 				end
% 				catdata=catdata';
% 				plot(obj.shiftIn10Hz_vec/10,catdata,Color=[1,1,1]*0.7);
% 				hold on;
% 				plotMeanSE(catdata,1,x=obj.shiftIn10Hz_vec/10,Color="blue",ShadeAlpha=0.5);
% 				titletext=sprintf("%s, # of cells = %g (%.2g%s)",titletext,size(catdata,2),size(catdata,2)/size(totalLog,1)*100,"%");
% 				title(titletext);
% 			end
% 		end
% 		%%
% 		function plotRidgeChangeAcrossShiftAllMice_ShiftTimePointsOfCellData(obj,ridgeparam2plot,options)
% 			arguments
% 				obj
% 				ridgeparam2plot=[]
% 				options.SkewThresh=[]
% 				options.CorrThresh=[]
% 			end
% 			flds=string(fieldnames(obj.ridgeResults_setting_plane_CellRidgeShift));
% 			alldata=struct2cell(obj.ridgeResults_setting_plane_CellRidgeShift);
% 			obj.plotTile_allPlane(@plotfunc,alldata,regexprep(flds,"_"," "));
% 			ChangeFigureSettings(XLabel="Time shift (sec)",Ylabel="Correlation",XLabelSize=20,YlabelSize=20, ...
% 				TitleSize=20,MasterTitle=sprintf("Correlation coefficient over shifting data, Skew>%g, maxCorr>%g, All mice combined", ...
% 				options.SkewThresh,options.CorrThresh),MasterTitleSize=20, ...
% 				Legend=[obj.info.NoAltplane.nameFull_p;"Std Error";"Mean"],LegendLoc="east",LegendSize=10);
% 
% 			function plotfunc(cdata,titletext)
% 				Skewlog=generateSkewList(obj,options.SkewThresh);
% 				catdata=[];
% 				totalLog=[];
% 				for i=1:length(cdata)
% 					currdata=cdata{i};
% 					if isempty(ridgeparam2plot)
% 						data=max(currdata,[],2);
% 					else
% 						IND=createIND(currdata,2,find(obj.ridgeparams==ridgeparam2plot));
% 						data=currdata(IND{:});
% 					end
% 					data=DeleteSingleDim(data,2);
% 					sLog=Skewlog{i};
% 					if ~isempty(options.CorrThresh)
% 						temp=max(data,[],2);
% 						cLog=temp>=options.CorrThresh;
% 					else
% 						cLog=ones(size(data,1),1,"logical");
% 					end
% 					Log=cLog(:).*sLog(:);
% 					data=data(logical(Log),:);
% 					data=mean(data,1);
% 					catdata=cat(1,catdata,data);
% 					totalLog=cat(1,totalLog,Log);
% 				end
% 				catdata=catdata';
% 				h=plot(obj.shiftIn10Hz_vec/10,catdata,Color=[1,1,1]*0.7);
% 				set(h,{'Color'},num2cell(parula(size(catdata,2)),2));
% 				hold on;
% 				plotMeanSE(catdata,1,x=obj.shiftIn10Hz_vec/10,Color="blue",ShadeAlpha=0.3);
% 				titletext=sprintf("%s, # of cells = %g (%.2g%s)",titletext,sum(totalLog),sum(totalLog)/size(totalLog,1)*100,"%");
% 				title(titletext);
% 			end
% 		end
% 		%%
% 		function plotRidgeChangeAcrossShiftAllMice_BrainRegions_ShiftTimePoints(obj,ridgeparam2plot,options)
% 			arguments
% 				obj
% 				ridgeparam2plot=[]
% 				options.SkewThresh=[]
% 				options.CorrThresh=[]
% 			end
% 			flds=string(fieldnames(obj.ridgeResults_setting_plane_CellRidgeShift));
% 			alldata=struct2cell(obj.ridgeResults_setting_plane_CellRidgeShift);
% 			ROIlog=generateROIlog(obj);
% 			obj.plotTile_allPlane(@plotfunc,alldata,regexprep(flds,"_"," "),1:length(alldata),TileHandle=gcf().Children,Exec_nexttile=false);
% 			ChangeFigureSettings(XLabel="Time shift (sec)",Ylabel="Correlation",XLabelSize=13,YlabelSize=13, ...
% 				TitleSize=20,MasterTitle=sprintf("Correlation coefficient over shifting data, Skew>%g, maxCorr>%g, All mice combined", ...
% 				options.SkewThresh,options.CorrThresh),MasterTitleSize=20, ...
% 				Legend=[obj.mouseName_MouseIND,"Std Error","Mean"],LegendLoc="east",LegendSize=10);
% 
% 			function plotfunc(cdata,titletext,INDex)
% 				Skewlog=generateSkewList(obj,options.SkewThresh);
% 				data_reg=[];
% 				for j=1:length(obj.contour_LeftRegion)	% for all regions
% 					catdata=[];
% 					totalLog=[];
% 					for i=1:length(cdata)	% planes
% 						currdata=cdata{i};
% 						if isempty(ridgeparam2plot)
% 							data=max(currdata,[],2);
% 						else
% 							IND=createIND(currdata,2,ridgeparam2plot);
% 							data=data(IND{:});
% 						end
% 						data=DeleteSingleDim(data,2);
% 
% 						% cell brain regions
% 						regLog=ROIlog{i}{j};
% 
% 
% 						sLog=Skewlog{i};
% 						if ~isempty(options.CorrThresh)
% 							temp=max(data,[],2);
% 							cLog=temp>=options.CorrThresh;
% 						else
% 							cLog=ones(size(data,1),1,"logical");
% 						end
% 						Log=cLog(:).*sLog(:).*regLog;
% 						data=data(logical(Log),:);
% 						data=mean(data,1);
% 						catdata=cat(1,catdata,data);
% 						totalLog=cat(1,totalLog,Log);
% 					end
% 					catdata=catdata';	% shift x plane	
% 
% 					% combine hemispheres
% 					combdata=[];
% 					for i=1:length(obj.plane_MouseIND)
% 						ind=obj.plane_MouseIND{i};
% 						temp=[];
% 						for k=1:length(ind)
% 							temp=cat(2,temp,catdata(:,ind(k)));
% 						end
% 						temp=mean(temp,2,"omitnan");
% 						combdata(:,i)=temp;
% 					end
% 
% 					data_reg{j}=combdata;	% shift x mouse
% 				end
% 
% 				fig=gcf;
% 				t=fig.Children;
% 				T=tiledlayout(t,'flow');
% 				T.Layout.Tile=INDex;
% 				T.Title.String=titletext;
% 				titletext2=convertCharsToStrings(obj.RegionNames_type);
% 				obj.plotTile_allPlane(@plotfunc2,data_reg,titletext2,TileHandle=T);
% 			end
% 			function plotfunc2(catdata,titletext2)
% 				h=plot(obj.shiftIn10Hz_vec/10,catdata,Color=[1,1,1]*0.7);
% 				set(h,{'Color'},num2cell(parula(size(catdata,2)),2));
% 				hold on;
% 				Log=all(isnan(catdata),1);
% 				catadata_nanomit=catdata(:,~Log);
% 				plotMeanSE(catadata_nanomit,1,x=obj.shiftIn10Hz_vec/10,Color="blue",ShadeAlpha=0.3);
% 				title(titletext2)
% 			end
% 		end
% 		%%
% 		function Skewlog=generateSkewList(obj,SkewnessThresh)
% 			Skewlog=[];
% 			for i=1:length(obj.info.NoAltplane.nameFull_p)
% 				csp=copy(obj);
% 				csp.plane=i;
% 				csp.load_specified(csp.plane,"Skew");
% 				% filter out cells with low skewness
% 				if ~isempty(SkewnessThresh)
% 					Skewlog{i}=csp.Skew(:)>=SkewnessThresh;
% 				else
% 					Skewlog{i}=logical(ones(size(csp.Skew(:))));
% 				end
% 			end
% 		end
% 		function ROIlog=generateROIlog(obj)
% 			ROIlog=[];
% 			for i=1:length(obj.info.NoAltplane.nameFull_p)
% 				csp=copy(obj);
% 				csp.plane=i;
% 				csp.load_specified(csp.plane,"regROI");
% 				ROIlog{i}=csp.Reg_inRegLog_reg;
% 			end
% 		end
% 	end
% 
% 
% 
% 
% 	methods
% 		%%
% 		function plotCortexDistributionBothHemisphere_ridge_Subplot(obj,settings,options)
% 			arguments
% 				obj
% 				settings.ridgeSetting2plot		string
% 				settings.ridgeparam2plot		 =[]
% 				settings.skew_thresh			(1,1) double
% 				settings.corr_thresh			(1,1) double
% 				settings.mask_radius_in_um		(1,1) double
% 				settings.mask_shape				(1,1) string {mustBeMember(settings.mask_shape, ["Gaussian","Ellipse"])}
% 				settings.padding_in_mm			(1,1) double = 0
% 				settings.FontSize=[]
% 
% 				options.overlayColormap {mustBeNumeric} = jet(256)
% 				options.alphaVal (1,1) double {mustBeGreaterThanOrEqual(options.alphaVal,0), ...
% 					mustBeLessThanOrEqual(options.alphaVal,1)} = 1
% 				options.mask_thresh {mustBeNumeric} = []
% 				options.Title (1,:) char = ''
% 				options.imagescmax (1,:) double = 1
% 				options.Plot_mask4conv (1,1) logical=true
% 			end
% 			ridgeSetting2plot = settings.ridgeSetting2plot;
% 			ridgeparam2plot    = settings.ridgeparam2plot;
% 			skew_thresh        = settings.skew_thresh;
% 			corr_thresh        = settings.corr_thresh;
% 			mask_radius_in_um  = settings.mask_radius_in_um;
% 			mask_shape         = settings.mask_shape;
% 			padding_in_mm      = settings.padding_in_mm;
% 
% 			ridgeResults_cell=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot);
% 			ridgeparams=obj.ridgeparams;
% 			obj.compareLength_warningUponMismatch(ridgeResults_cell,obj.info.NoAltplane.nameFull_p);
% 			halfplanevec=1:length(ridgeResults_cell)/2;
% 			mask=CSP.generate_mask4conv_in_plotCellDensity(mask_radius_in_um,mask_shape,padding_in_mm);
% 			imagescmax=options.imagescmax;
% 			options.AtlasImage=obj.AtlasImage;
% 			ridgeResults_cell1=ridgeResults_cell(1:2:end);
% 			ridgeResults_cell2=ridgeResults_cell(2:2:end);
% 			nameMouse=obj.info.NoAltplane.nameMouse_p(1:2:end);
% 			obj.plotSubPlot_allPlane(@plotfunc,ridgeResults_cell1,ridgeResults_cell2,halfplanevec,nameMouse,linkaxes="xy");
% 			sgt=sgtitle(sprintf("Density of cells with high correlation for model: %s. Corr threshold: %.2g, " + ...
% 				"\nSkew threshold: %.2g, MaskRadius: %.3g \\mum, Cell count threshold: %g", ...
% 				regexprep(ridgeSetting2plot,"_"," "),corr_thresh,skew_thresh,mask_radius_in_um,options.mask_thresh));
% 			try; sgt.FontSize=settings.FontSize; end
% 
% 			function plotfunc(ridgeResults1,ridgeResults2,halfplane,titletext)
% 				obj.plane=halfplane*2-1;
% 				obj.load_specified(obj.plane,"Skew");
% 				obj.load_specified(obj.plane,"regROI");
% 				% filter out cells with low skewness
% 				Skewlog=obj.Skew(:)>=skew_thresh;
% 				ridgeResults=ridgeResults1(Skewlog,:);
% 				regROI=obj.regROI(Skewlog,:);
% 				% obtain cell log for correlation
% 				if isempty(ridgeparam2plot)
% 					[corrval,corrind]=max(ridgeResults,[],2);
% 					cellLog=corrval>=corr_thresh;
% 				else
% 					paramIND=find(ridgeparams==ridgeparam2plot,1);
% 					cellLog=ridgeResults(:,paramIND)>=corr_thresh;
% 				end
% 				name1=obj.info.NoAltplane.nameMouse_p(obj.plane);
% 				regROI1=regROI;
% 				cellLog1=cellLog;
% 
% 				obj.plane=halfplane*2;
% 				obj.load_specified(obj.plane,"Skew");
% 				obj.load_specified(obj.plane,"regROI");
% 				% filter out cells with low skewness
% 				Skewlog=obj.Skew(:)>=skew_thresh;
% 				ridgeResults=ridgeResults2(Skewlog,:);
% 				regROI=obj.regROI(Skewlog,:);
% 				% obtain cell log for correlation
% 				if isempty(ridgeparam2plot)
% 					[corrval,corrind]=max(ridgeResults,[],2);
% 					cellLog=corrval>=corr_thresh;
% 				else
% 					paramIND=find(ridgeparams==ridgeparam2plot,1);
% 					cellLog=ridgeResults(:,paramIND)>=corr_thresh;
% 				end
% 				name2=obj.info.NoAltplane.nameMouse_p(obj.plane);
% 				regROI2=regROI;
% 				cellLog2=cellLog;
% 
% 				if name1~=name2;error("mice names must match");end
% 
% 				regROI=[regROI1;regROI2];
% 				cellLog=[cellLog1;cellLog2];
% 
% 				try;options.imagescmax=imagescmax(halfplane);end
% 				opts=namedargs2cell(options);
% 				plotCellDensity(regROI, cellLog, mask, opts{:})
% 				title(sprintf("%s, cells with Skew>%.2g: %g,\ncells with Corr>=%2g: %g, (%.1f%s",titletext,skew_thresh,length(cellLog),corr_thresh,sum(cellLog),sum(cellLog)/length(cellLog)*100,"%)"));
% 			end
% 		end
% 		function plotCortexDistributionBothHemisphere_ridge_allPlanes(obj,settings,options)
% 			arguments
% 				obj
% 				settings.ridgeSetting2plot		string
% 				settings.ridgeparam2plot		(1,1) double
% 				settings.skew_thresh			(1,1) double
% 				settings.corr_thresh			(1,1) double
% 				settings.mask_radius_in_um		(1,1) double
% 				settings.mask_shape				(1,1) string {mustBeMember(settings.mask_shape, ["Gaussian","Ellipse"])}
% 				settings.padding_in_mm			(1,1) double = 0
% 				settings.FontSize=[]
% 
% 				options.overlayColormap {mustBeNumeric} = jet(256)
% 				options.alphaVal (1,1) double {mustBeGreaterThanOrEqual(options.alphaVal,0), ...
% 					mustBeLessThanOrEqual(options.alphaVal,1)} = 1
% 				options.mask_thresh {mustBeNumeric} = []
% 				options.Title (1,:) char = ''
% 				options.imagescmax (1,:) double = 1
% 				options.Plot_mask4conv (1,1) logical=true
% 			end
% 			ridgeSetting2plot = settings.ridgeSetting2plot;
% 			ridgeparam2plot    = settings.ridgeparam2plot;
% 			skew_thresh        = settings.skew_thresh;
% 			corr_thresh        = settings.corr_thresh;
% 			mask_radius_in_um  = settings.mask_radius_in_um;
% 			mask_shape         = settings.mask_shape;
% 			padding_in_mm      = settings.padding_in_mm;
% 
% 			ridgeparams=obj.ridgeparams;
% 			ridgeResults_cell=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot);
% 			obj.compareLength_warningUponMismatch(ridgeResults_cell,obj.info.NoAltplane.nameFull_p);
% 			mask=CSP.generate_mask4conv_in_plotCellDensity(mask_radius_in_um,mask_shape,padding_in_mm);
% 			allCellLog=[];
% 			allregROI=[];
% 			for i=1:length(ridgeResults_cell)
% 				ridgeResults=ridgeResults_cell{i};
% 				if isempty(ridgeparam2plot)
% 					[corrval,corrind]=max(ridgeResults,[],2);
% 					ridgeResults=corrval;
% 				else
% 					paramIND=find(ridgeparams==ridgeparam2plot,1);
% 					ridgeResults=ridgeResults(:,paramIND);
% 				end
% 				obj.plane=i;
% 				obj.load_specified(obj.plane,"Skew");
% 				obj.load_specified(obj.plane,"regROI");
% 				Skewlog=obj.Skew(:)>=skew_thresh;
% 				ridgeResults=ridgeResults(Skewlog);
% 				regROI=obj.regROI(Skewlog,:);
% 				cellLog=ridgeResults>=corr_thresh;
% 				allCellLog=cat(1,allCellLog,cellLog);
% 				allregROI=cat(1,allregROI,regROI);
% 			end
% 			options.AtlasImage=obj.AtlasImage;
% 			opts=namedargs2cell(options);
% 			plotCellDensity(allregROI, allCellLog, mask, opts{:});
% 			htitle=title(sprintf("Density of cells with high correlation for model: %s. Corr threshold: %.2g, \nSkew threshold: %.2g, MaskRadius: %.3g \\mum, Cell count threshold: %g", ...
% 				regexprep(ridgeSetting2plot,"_"," "),corr_thresh,skew_thresh,mask_radius_in_um,options.mask_thresh),FontSize=settings.FontSize);
% 			% title(sprintf("%s, cells with Skew>%.2g: %g,\ncells with Corr>=%2g: %g, (%.1f%s",titletext,skew_thresh,length(cellLog),corr_thresh,sum(cellLog),sum(cellLog)/length(cellLog)*100,"%)"));
% 			% set(htitle, 'Clipping', 'off')
% 			ax = gca;
% 			pos = ax.Position;
% 			% Decrease the vertical position and height to leave more room on top.
% 			ax.Position = [pos(1) pos(2)-0.05 pos(3) pos(4)*0.95];
% 
% 		end
% 
% 		function plotCortexDistributionComparisonBothHemisphere_ridge_allPlanes(obj,settings,options)
% 			arguments
% 				obj
% 				settings.ridgeSetting_target
% 				settings.ridgeSetting_base4comparison
% 				settings.func_handle_outputMustBeLogical
% 				settings.ridgeparam2plot		=[]
% 				settings.skew_thresh			(1,1) double
% 				settings.corr_thresh			(1,1) double
% 				settings.mask_radius_in_um		(1,1) double
% 				settings.mask_shape				(1,1) string {mustBeMember(settings.mask_shape, ["Gaussian","Ellipse"])}
% 				settings.padding_in_mm			(1,1) double = 0
% 
% 				options.overlayColormap {mustBeNumeric} = jet(256)
% 				options.alphaVal (1,1) double {mustBeGreaterThanOrEqual(options.alphaVal,0), ...
% 					mustBeLessThanOrEqual(options.alphaVal,1)} = 1
% 				options.mask_thresh {mustBeNumeric} = []
% 				options.Title (1,:) char = ''
% 				options.imagescmax (1,:) double = 1
% 				options.Plot_mask4conv (1,1) logical=true
% 			end
% 			ridgeSetting_target=settings.ridgeSetting_target;
% 			ridgeSetting_base=settings.ridgeSetting_base4comparison;
% 			func_handle=settings.func_handle_outputMustBeLogical;
% 			ridgeparam2plot    = settings.ridgeparam2plot;
% 			skew_thresh        = settings.skew_thresh;
% 			corr_thresh        = settings.corr_thresh;
% 			mask_radius_in_um  = settings.mask_radius_in_um;
% 			mask_shape         = settings.mask_shape;
% 			padding_in_mm      = settings.padding_in_mm;
% 			ridgeparams=obj.ridgeparams;
% 			mask=CSP.generate_mask4conv_in_plotCellDensity(mask_radius_in_um,mask_shape,padding_in_mm);
% 
% 			[allRidgeResults_target,allregROI_target,allCellLog_target]=concatenateAllData(ridgeSetting_target);
% 			[allRidgeResults_base,allregROI_base,allCellLog_base]=concatenateAllData(ridgeSetting_base);
% 			allcustomCellLog=func_handle(allRidgeResults_target,allRidgeResults_base);
% 			if ~isequal(allregROI_target,allregROI_base);error("ROI does not match");end
% 
% 			options.AtlasImage=obj.AtlasImage;
% 			opts=namedargs2cell(options);
% 			plotCellDensity(allregROI_target, allcustomCellLog, mask, opts{:});
% 
% 			strfunc=generateFuncHandleString(func_handle);
% 			title(sprintf("Density of cells satisfying the condition (%s). \nCorr threshold: %.2g, Skew threshold: %.2g, MaskRadius: %.3g mm, Cell count threshold: %g", ...
% 				strfunc,corr_thresh,skew_thresh,mask_radius_in_um,options.mask_thresh),FontSize=15);
% 			% title(sprintf("%s, cells with Skew>%.2g: %g,\ncells with Corr>=%2g: %g, (%.1f%s",titletext,skew_thresh,length(cellLog),corr_thresh,sum(cellLog),sum(cellLog)/length(cellLog)*100,"%)"));
% 			ax = gca;
% 			pos = ax.Position;
% 			% Decrease the vertical position and height to leave more room on top.
% 			ax.Position = [pos(1) pos(2)-0.05 pos(3) pos(4)*0.95];
% 
% 			function strfunc=generateFuncHandleString(func_handle)
% 				strfunc=string(str2sym(func2str(func_handle)));
% 				strfunc=regexprep(strfunc,"target",sprintf("[%s]",ridgeSetting_target));
% 				strfunc=regexprep(strfunc,"base",sprintf("[%s]",ridgeSetting_base));
% 				strfunc=regexprep(strfunc,"corr_thresh",sprintf("%g",corr_thresh));
% 				strfunc=regexprep(strfunc,"_"," ");
% 			end
% 			function [allRidgeResults,allregROI,allCellLog]=concatenateAllData(ridgeSetting2plot)
% 				ridgeResults_cell=obj.ridgeResults_setting_plane_CellRidge.(ridgeSetting2plot);
% 				obj.compareLength_warningUponMismatch(ridgeResults_cell,obj.info.NoAltplane.nameFull_p);
% 				mask=CSP.generate_mask4conv_in_plotCellDensity(mask_radius_in_um,mask_shape,padding_in_mm);
% 				allRidgeResults=[];
% 				allCellLog=[];
% 				allregROI=[];
% 				for i=1:length(ridgeResults_cell)
% 					ridgeResults=ridgeResults_cell{i};
% 					if isempty(ridgeparam2plot)
% 						[corrval,corrind]=max(ridgeResults,[],2);
% 						ridgeResults=corrval;
% 					else
% 						paramIND=find(ridgeparams==ridgeparam2plot,1);
% 						ridgeResults=ridgeResults(:,paramIND);
% 					end
% 					obj.plane=i;
% 					obj.load_specified(obj.plane,"Skew");
% 					obj.load_specified(obj.plane,"regROI");
% 					Skewlog=obj.Skew(:)>=skew_thresh;
% 					ridgeResults=ridgeResults(Skewlog);
% 					regROI=obj.regROI(Skewlog,:);
% 					cellLog=ridgeResults>=corr_thresh;
% 					allRidgeResults=cat(1,allRidgeResults,ridgeResults);
% 					allCellLog=cat(1,allCellLog,cellLog);
% 					allregROI=cat(1,allregROI,regROI);
% 				end
% 			end
% 		end
% 
% 
% 		function [regROI,cellLog]=PreprocessData2plotCortexDistribution_ridge(obj,ridgeResults,corr_thresh,skewthresh,options)
% 			% regROI,cellLog are both vertical vectors (filtered by skew)
% 			arguments
% 				obj
% 				ridgeResults
% 				corr_thresh
% 				skewthresh
% 				options.ridgeparams
% 				options.ridgeparam2use
% 			end
% 			obj.load_specified(obj.plane,"Skew");
% 			obj.load_specified(obj.plane,"regROI");
% 			Skewlog=obj.Skew(:)>=skewthresh;
% 			ridgeResults=ridgeResults(Skewlog,:);
% 			regROI=obj.regROI(Skewlog,:);
% 
% 			[corrval,corrind]=max(ridgeResults,[],2);
% 			cellLog=corrval>=corr_thresh;
% 		end
% 		function [regROI,cellLog]=PreprocessData2plotCortexDistribution_ridge_outdated(obj,ridgeResults,corr_thresh,skewthresh,options)
% 			% regROI,cellLog are both vertical vectors (filtered by skew)
% 			arguments
% 				obj
% 				ridgeResults
% 				corr_thresh
% 				skewthresh
% 				options.ridgeparams
% 				options.ridgeparam2use
% 			end
% 			obj.load_specified(obj.plane,"Skew");
% 			obj.load_specified(obj.plane,"regROI");
% 			Skewlog=obj.Skew(:)>=skewthresh;
% 			ridgeResults=ridgeResults(Skewlog,:);
% 			regROI=obj.regROI(Skewlog,:);
% 
% 			[corrval,corrind]=max(ridgeResults,[],2);
% 			cellLog=corrval>=corr_thresh;
% 		end
% 		% function plotCortexDistribution_ridge(obj,ridgeResults,corr_thresh,ridgeparams,ridgeparam2use,radius_in_um,option,padding,options)
% 		% 	arguments
% 			% 	obj
% 			% 	ridgeResults
% 			% 	corr_thresh
% 			% 	ridgeparams
% 			% 	ridgeparam2use
% 			% 	radius_in_um
% 			% 	option (1,1) string {mustBeMember(option, ["Gaussian","Ellipse"])}
% 			% 	padding.padding_in_mm=0
% 		% 
% 			% 	options.overlayColormap {mustBeNumeric} = jet(256)
% 			% 	options.alphaVal (1,1) double {mustBeGreaterThanOrEqual(options.alphaVal,0), ...
% 				% 	mustBeLessThanOrEqual(options.alphaVal,1)} = 1
% 			% 	options.mask_thresh {mustBeNumeric} = []
% 			% 	options.Title (1,:) char = ''
% 			% 	options.imagescmax (1,1) double = 1
% 			% 	options.Plot_mask4conv (1,1) logical=true
% 		% 	end
% 		% 	[corrval,corrind]=max(ridgeResults,[],2);
% 		% 	cellLog=corrval>=corr_thresh;
% 		% 	% obtain cell regROI
% 		% 	obj.load_specified(obj.plane,"regROI");
% 		% 	padding=namedargs2cell(padding);
% 		% 	mask=CSP.generate_mask4conv_in_plotCellDensity(radius_in_um,option,padding{:});
% 		% 
% 		% 	options.AtlasImage=obj.AtlasImage;
% 		% 	options=namedargs2cell(options);
% 		% 	CSP.plotCellDensityUsingMaskConvolution(obj.regROI,cellLog,mask,options{:});
% 		% end
% 	end















































%% combine normal SP and DLC data
properties
	SP_and_DLC_combined
	SP_and_DLC_combined_names
end
methods
	function combine_SP_and_DLC(obj,options)
		arguments
			obj
			options.SP (1,:) string {mustBeMember(options.SP,["theta","phi","tau","X","Y","Z",""])}
			options.DLC (1,:) string {mustBeMember(options.DLC,["tXYZbp","tDXYZbp","tD2XYZbp", ...
				"tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp", ...
				"tDXYZbp_magnitude","tDDOTbp_magnitude","tDANGbp_magnitude", ...
				"tD2XYZbp_magnitude","tD2DOTbp_magnitude","tD2ANGbp_magnitude", ...
				"tDXYZbp_magnitude_mean", ...
				"eigen_DLC",""])}
		end
		XYZ=["tXYZbp","tDXYZbp","tD2XYZbp"];
		Ang=["tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"];
		SPnamelist=string(options.SP);
		DLCnamelist=string(options.DLC);

		SPdata=[];
		SPnames=[];
		if options.SP~=""
			SPind=obj.ind_SPparam(SPnamelist);
			ind=makeIND(obj.tpk,2,SPind);
			SPdata=obj.tpk(ind{:});
			SPnames=SPnamelist;
		end

		DLCdata=[];
		DLCnames=[];
		if options.DLC~=""
			obj.make_DLC_speed;
			obj.make_DLC_accel;
			for i=1:length(DLCnamelist)
				currstr=DLCnamelist(i);
				name_of_current_category=[];
				if isfield(obj.DLC,currstr)
					numDLC=obj.DLC_convert2numeric(obj.DLC);
					currdata=numDLC.(currstr);
					if any(containsSubstring(XYZ,currstr))
						temp=arrayfun(@(a)currstr+"_"+a+["_x","_y","_z"],obj.rec3d.bp_ind,'UniformOutput',false);
						name_of_current_category=cat(2,temp{:});...1 x 42 string array
					elseif any(containsSubstring(Ang,currstr))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif isfield(obj.DLC_speed,currstr)
					currdata=obj.DLC_speed.(currstr);
					if any(containsSubstring(currstr,"XYZ"))
						name_of_current_category=currstr+"_"+obj.rec3d.bp_ind;...1 x 14 string array
					elseif any(containsSubstring(currstr,"ANG"))||any(containsSubstring(currstr,"DOT"))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif isfield(obj.DLC_accel,currstr)
					currdata=obj.DLC_accel.(currstr);
					if any(containsSubstring(currstr,"XYZ"))
						name_of_current_category=currstr+"_"+obj.rec3d.bp_ind;...1 x 14 string array
					elseif any(containsSubstring(currstr,"ANG"))||any(containsSubstring(currstr,"DOT"))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif currstr=="tDXYZbp_magnitude_mean"
					data=obj.DLC_speed.tDXYZbp_magnitude;
					currdata=mean(data,2);
					name_of_current_category="tDXYZbp_magnitude_mean";
				elseif currstr=="eigen_DLC"
					if isempty(obj.eigen_DLC)
						error("solve for eigenposture first")
					end
					currdata=obj.eigen_DLC;
					name_of_current_category="eigen_posture_"+string(1:size(currdata,2));
				else
					error("Specified param2discretize is not recognized")
				end
				DLCdata=cat(2,DLCdata,currdata);
				DLCnames=[DLCnames,name_of_current_category];
				if size(DLCdata,2)~=size(DLCnames,2)
					error("Size of DLCdata and DLCnames must match");
				end
			end
		end
		obj.SP_and_DLC_combined=[SPdata,DLCdata];
		obj.SP_and_DLC_combined_names=[SPnames,DLCnames];
		if size(obj.SP_and_DLC_combined,2)~=size(obj.SP_and_DLC_combined_names,2)
			error("Size of SP_and_DLC_combined and SP_and_DLC_combined_names must match");
		end
	end
end

%% combine Discretized SP and DLC data
properties
	Discretized_SP_and_DLC
	Discretized_SP_and_DLC_names
	Discretized_SP_and_DLC_undiscretizer
	Discretized_SP_and_DLC_undiscretized_names
	Discretized_SP_and_DLC_edges
end
methods
	function combine_Discretized_SP_and_DLC(obj,options)
		arguments
			obj
			options.Data2combine (1,1) string {mustBeMember(options.Data2combine, ["Both","DLC only","SP only"])}
			% options.CombineRawValues (1,1) logical =true
		end
		% params used
		%{
		SP_Discretized_cellVector
		DLC_Discretized_cellVector
		SP_Discretized_names			% cell: {["theta_phi"]}    {["theta"]}
		DLC_Discretized_names			% cell: {1×42 string}    {1×16  string}
		SP_Discretized_params2encode	% string array: "theta"    "phi"    "theta"
		%}
		SPdata=obj.SP_Discretized_cellVector;
		SPdata=cell2mat(SPdata(:)');
		DLCdata=obj.DLC_Discretized_cellVector;
		DLCdata=cell2mat([DLCdata{:}]);
		% categorize names
		SPdata_names=obj.SP_Discretized_names;
		SPdata_names=cat(2,SPdata_names{:});
		DLCdata_names=obj.DLC_Discretized_names;
		temp=[DLCdata_names{:}];
		% DLCdata_names=cat(2,temp{:});
		DLCdata_names=temp;
		if size(SPdata,2)~=size(SPdata_names,2);error("The size of array and names do not match");end
		if size(DLCdata,2)~=size(DLCdata_names,2);error("The size of array and names do not match");end
		
		% categorize Undiscretizer
		SP_Undiscretizer=obj.SP_Discretized_Undiscretizer(:)';
		DLC_Undiscretizer=[obj.DLC_Discretized_Undiscretizer{:}];

		% categorize undiscretized names
		SP_UndiscretizedNames=[obj.SP_Discretized_Undiscretized_names{:}];
		DLC_UndiscretizedNames=[obj.DLC_Discretized_Undiscretized_names{:}];

		% categorize edges
		SP_edges=obj.SP_Discretized_edges(:)';
		DLC_edges=[obj.DLC_Discretized_edges_cellVector{:}];

		% XYZ=["tXYZbp","tDXYZbp","tD2XYZbp"];
		% Ang=["tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"];		
		% if options.CombineRawValues
		% 	% SP data
		% 	currstr=obj.SP_Discretized_params2encode;
		% 	currstr=unique(currstr);
		% 	IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam(currstr);
		% 	SPdata=cat(2,SPdata,obj.tpk(IND{:}));
		% 	% concatenate names
		% 	currname=currstr;
		% 	SPdata_names=cat(2,SPdata_names,currname);
		% 	if size(SPdata,2)~=size(SPdata_names,2);error("The size of array and names do not match");end
		% 
		% 	% DLC data
		% 	switch obj.DLC_inv_or_normal
			% 	case "normalPosture"
				% 	currdata=obj.DLC_convert(obj.DLC,"AllNumeric");
			% 	case "PostureInversedBySP"
				% 	currdata=obj.DLC_convert(obj.invDLC,"AllNumeric");
		% 	end
		% 	DLC_normal=currdata;
		% 	currstr=obj.DLC_Discretized_params2encode;
		% 	for k=currstr(:)'
			% 	name_of_current_category=[];
			% 	if isfield(DLC_normal,k)
				% 	DLCdata=cat(2,DLCdata,DLC_normal.(k));
				% 	if any(containsSubstring(XYZ,k))
					% 	temp=arrayfun(@(a)k+"_"+a+["_x","_y","_z"],obj.rec3d.bp_ind,'UniformOutput',false);
					% 	name_of_current_category=cat(2,temp{:});...1 x 42 string array
				% 	elseif any(containsSubstring(Ang,k))
					% 	name_of_current_category=k+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
				% 	else
					% 	error("No match")
				% 	end
			% 	elseif isfield(obj.DLC_speed,k)
				% 	temp=obj.DLC_speed;
				% 	DLCdata=cat(2,DLCdata,temp.(k));
				% 	if any(containsSubstring(k,"XYZ"))
					% 	name_of_current_category=k+"_"+obj.rec3d.bp_ind;...1 x 14 string array
				% 	elseif any(containsSubstring(k,"ANG"))||any(containsSubstring(k,"DOT"))
					% 	name_of_current_category=k+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
				% 	else
					% 	error("No match")
				% 	end
			% 	elseif isfield(obj.DLC_accel,k)
				% 	temp=obj.DLC_accel;
				% 	DLCdata=cat(2,DLCdata,temp.(k));
				% 	if any(containsSubstring(k,"XYZ"))
					% 	name_of_current_category=k+"_"+obj.rec3d.bp_ind;...1 x 14 string array
				% 	elseif any(containsSubstring(k,"ANG"))||any(containsSubstring(k,"DOT"))
					% 	name_of_current_category=k+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
				% 	else
					% 	error("No match")
				% 	end
			% 	else
				% 	error("No match")
			% 	end
			% 	name_of_current_category=name_of_current_category+"_raw";
			% 	DLCdata_names=cat(2,DLCdata_names,name_of_current_category);
			% 	if size(DLCdata,2)~=size(DLCdata_names,2);error("The size of array and names do not match");end
		% 	end
		% end
		if size(SPdata,2)~=size(SPdata_names,2);error("The size of array and names do not match");end
		if size(DLCdata,2)~=size(DLCdata_names,2);error("The size of array and names do not match");end
		switch options.Data2combine
			case "Both"
				all_Pred=[SPdata,DLCdata];
				all_names=[SPdata_names,DLCdata_names];
				undiscretizer=[SP_Undiscretizer,DLC_Undiscretizer];
				undiscretized_names=[SP_UndiscretizedNames,DLC_UndiscretizedNames];
				edges=[SP_edges,DLC_edges];
			case "DLC only"
				all_Pred=[DLCdata];
				all_names=[DLCdata_names];
				undiscretizer=[DLC_Undiscretizer];
				undiscretized_names=[DLC_UndiscretizedNames];
				edges=[DLC_edges];
			case "SP only"
				all_Pred=[SPdata];
				all_names=[SPdata_names];
				undiscretizer=[SP_Undiscretizer];
				undiscretized_names=[SP_UndiscretizedNames];
				edges=[SP_edges];	
		end
		obj.Discretized_SP_and_DLC=all_Pred;
		obj.Discretized_SP_and_DLC_names=all_names;
		obj.Discretized_SP_and_DLC_undiscretizer=undiscretizer;
		obj.Discretized_SP_and_DLC_undiscretized_names=undiscretized_names;
		obj.Discretized_SP_and_DLC_edges=edges;
	end
end




%% split data to x validation groups for training and testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
	xval_Dataset
end
methods
	function split_Data_to_xvalidation_groups(obj,name,array,options)
		arguments
			obj
		end
		arguments(Repeating)
			name (1,1) string
			array 
		end
		arguments
			options.Predictor	% time x param x trial x ...
			options.Responser	% time x param x trial x ...
			options.groupnum
			options.dim2split=3
		end
		Trainer=[];
		Tester=[];
		obj.groupnum=options.groupnum;
		dim2split=options.dim2split;
		% check dim
		for dim=1:ndims(options.Predictor)
			dimlist=[];
			if dim==options.dim2split
				dimlist(end+1)=size(options.Predictor,dim2split);
				dimlist(end+1)=size(options.Responser,dim2split);
				if ~isempty(name)
					for i=1:length(array)
						dimlist(end+1)=size(array{i},dim2split);
					end
				end
				if numel(unique(dimlist))~=1
					error("the size along dim2split has to be the same")
				end
			elseif dim~=2
				dimlist(end+1)=size(options.Predictor,dim2split);
				dimlist(end+1)=size(options.Responser,dim2split);
				if ~isempty(name)
					for i=1:length(array)
						dimlist(end+1)=size(array{i},dim2split);
					end
				end
				if numel(unique(dimlist))~=1
					warning(sprintf("the size along dim %d does not match",dim))
				end
			end
		end
		[Trainer.Predictor,Tester.Predictor]=splitArrayToGroups(options.Predictor, options.groupnum, dim2split);
		[Trainer.Responser,Tester.Responser]=splitArrayToGroups(options.Responser, options.groupnum, dim2split);
		for i=1:numel(array)
			[Trainer.(name{i}),Tester.(name{i})]=splitArrayToGroups(array{i},options.groupnum,dim2split);
		end
		obj.xval_Dataset=[];
		obj.xval_Dataset.Trainer=Trainer;
		obj.xval_Dataset.Tester=Tester;
		obj.xval_Dataset=StructWithCell2CellWithStruct_multiLayers(obj.xval_Dataset);
	end
	function split_Data_to_xvalidation_groups_continuous(obj,options)
		arguments
			obj
			options.Predictor	% time x param
			options.Responser	% time x param
			options.Others	=[]	% cell containing time x param
			options.groupnum
			options.TimeChunkSec (1,1) =5	% number of time points per chunk
		end
		Trainer=[];
		Tester=[];
		obj.groupnum=options.groupnum;
		pred=options.Predictor;
		resp=options.Responser;
		gnum=options.groupnum;
		chunk=options.TimeChunkSec*obj.SamplingRate_Hz;
		numchunk=floor(size(pred,1)./chunk);
		pred=pred(1:numchunk*chunk,:);
		resp=resp(1:numchunk*chunk,:);
		pred_tpk=reshape(pred',size(pred,2),chunk,numchunk);	pred_tpk=permute(pred_tpk,[2,1,3]);
		resp_tpk=reshape(resp',size(resp,2),chunk,numchunk);	resp_tpk=permute(resp_tpk,[2,1,3]);

		[Trainer.Predictor,Tester.Predictor]=splitArrayToGroups(pred_tpk, gnum, 3);
		[Trainer.Responser,Tester.Responser]=splitArrayToGroups(resp_tpk, gnum, 3);
		obj.xval_Dataset=[];
		obj.xval_Dataset.Trainer=Trainer;
		obj.xval_Dataset.Tester=Tester;
		obj.xval_Dataset=StructWithCell2CellWithStruct_multiLayers(obj.xval_Dataset);
	end
	function name=make_xval_name(obj,propStrArr,name2append)
		props=cellstr(propStrArr);
		props=cellfun(@string,props,UniformOutput=false);
		name=props;
		name{end}=name{end}+name2append;
	end
	% function set_xval_prop(obj,props,val)
	% 	props=string(props);
	% 	props=cellstr(props);
	% 	props=cellfun(@string,props,UniformOutput=false);
	% end
end
methods (Static)
	function newmat=reshape_mat_to_3rd_Dim(mat,options)
		arguments
			mat
			options.groupnum
			options.dim
		end
		groupnum=options.groupnum;
		dim=options.dim;
		if dim==2
			mat=mat';
		end
		timelength=floor(size(mat,1)/groupnum);
		mat=mat(1:timelength*groupnum,:);
		newmat=permute(reshape(mat',size(mat,2),timelength,groupnum),[2,1,3]);
		if dim==2
			newmat=pagetranspose(newmat);
		end
	end
end

%% dimensionality reduction %%%%%%%%%%%
methods
	function pca_xval_Dataset(obj,options)
		arguments
			obj
			options.Prop2ApplyPCA (1,:) string
			options.Prop2_Multiply_PCA_Axes (1,:) string =[]
			options.NumDim2keep=[]
		end
		Axes_name=makeName(options.Prop2ApplyPCA,"_pcaAxes");
		pca_data_name=makeName(options.Prop2ApplyPCA,"_pca");
		name1=makeName(options.Prop2ApplyPCA,"");
		for i=1:numel(obj.xval_Dataset)
			data=getfield(obj.xval_Dataset{i},name1{:});
			f=@(a)multiBlockCat(a,[3:ndims(a)],1);
			TKDc=f(data);
			Axes=pca(TKDc);
			numdim=size(Axes,2);
			if ~isempty(options.NumDim2keep)
				numdim=options.NumDim2keep;
				numdim=min(numdim,size(Axes,2));
			end
			IND=makeIND(Axes,2,1:numdim);
			Axes=Axes(IND{:});
			obj.xval_Dataset{i}=setfield(obj.xval_Dataset{i},Axes_name{:},Axes);
			obj.xval_Dataset{i}=setfield(obj.xval_Dataset{i},pca_data_name{:},pagemtimes(data,Axes));
			if ~isempty(options.Prop2_Multiply_PCA_Axes)
				pca_data_name2=makeName(options.Prop2_Multiply_PCA_Axes,"_pca");
				name2=makeName(options.Prop2_Multiply_PCA_Axes,"");
				data2=getfield(obj.xval_Dataset{i},name2{:});
				obj.xval_Dataset{i}=setfield(obj.xval_Dataset{i},pca_data_name2{:},pagemtimes(data2,Axes));
			end
		end
		function name=makeName(propStrArr,name2append)
			props=cellstr(propStrArr);
			props=cellfun(@string,props,UniformOutput=false);
			name=props;
			name{end}=name{end}+name2append;
		end
	end
	%%
	function cPCA_xval_Dataset(obj,options)
		arguments
			obj
			options.weight
			options.Prop2ApplycPCA (1,:) string
			options.Prop2_Multiply_cPCA_Axes (1,:) string =[]
			options.NumDim2keep=[]
		end
		Axes_name=makeName(options.Prop2ApplycPCA,"_cpcaAxes");
		cPCA_data_name=makeName(options.Prop2ApplycPCA,"_cpca");
		Prop2ApplycPCA=makeName(options.Prop2ApplycPCA,"");
		for i=1:numel(obj.xval_Dataset)
			data=getfield(obj.xval_Dataset{i},Prop2ApplycPCA{:});
			f=@(a)multiBlockCat(a,[4:ndims(a)],1);
			TDck=f(data);
			[Axes,Score,A,B,ScoreA,ScoreB,val]=nFDA(TDck,'nFDA',options.weight);
			numdim=size(Axes,2);
			if ~isempty(options.NumDim2keep)
				numdim=options.NumDim2keep;
				numdim=min(numdim,size(Axes,2));
			end
			IND=makeIND(Axes,2,1:numdim);
			Axes=Axes(IND{:});
			obj.xval_Dataset{i}=setfield(obj.xval_Dataset{i},Axes_name{:},Axes);
			obj.xval_Dataset{i}=setfield(obj.xval_Dataset{i},cPCA_data_name{:},pagemtimes(data,Axes));
			if ~isempty(options.Prop2_Multiply_cPCA_Axes)
				pca_data_name2=makeName(options.Prop2_Multiply_cPCA_Axes,"_cpca");
				Prop2_Multiply_cPCA_Axes=makeName(options.Prop2_Multiply_cPCA_Axes,"");
				data2=getfield(obj.xval_Dataset{i},Prop2_Multiply_cPCA_Axes{:});
				obj.xval_Dataset{i}=setfield(obj.xval_Dataset{i},pca_data_name2{:},pagemtimes(data2,Axes));
			end
		end
		function name=makeName(propStrArr,name2append)
			props=cellstr(propStrArr);
			props=cellfun(@string,props,UniformOutput=false);
			name=props;
			name{end}=name{end}+name2append;
		end
	end
	
end
%% CCA
methods (Static)
	function cca_xval_Dataset(csp1,csp2,options)
		arguments
			csp1 CSP
			csp2 CSP
			options.Prop2ApplyCCA (1,:) string
			options.Prop2_Multiply_CCA_Axes (1,:) string=[]
		end
		prop=makeName(options.Prop2ApplyCCA,"");
		ccaprop=makeName(options.Prop2ApplyCCA,"_cca");
		axprop=makeName(options.Prop2ApplyCCA,"_ccaAxes");
		corrprop=makeName(options.Prop2ApplyCCA,"_corr");

		for i=1:numel(csp1.xval_Dataset)
			data1=getfield(csp1.xval_Dataset{i},prop{:});
			data2=getfield(csp2.xval_Dataset{i},prop{:});
			f=@(a)multiBlockCat(a,[3:ndims(a)],1);
			TKDc1=f(data1);
			TKDc2=f(data2);
			[W1,W2,r,score1,score2]=canoncorr(TKDc1,TKDc2);
			csp1.xval_Dataset{i}=setfield(csp1.xval_Dataset{i},axprop{:},W1);
			csp2.xval_Dataset{i}=setfield(csp2.xval_Dataset{i},axprop{:},W2);
			csp1.xval_Dataset{i}=setfield(csp1.xval_Dataset{i},ccaprop{:},pagemtimes(data1,W1));
			csp2.xval_Dataset{i}=setfield(csp2.xval_Dataset{i},ccaprop{:},pagemtimes(data2,W2));
			csp1.xval_Dataset{i}=setfield(csp1.xval_Dataset{i},corrprop{:},r);
			csp2.xval_Dataset{i}=setfield(csp2.xval_Dataset{i},corrprop{:},r);
			if ~isempty(options.Prop2_Multiply_CCA_Axes)
				testccaprop=makeName(options.Prop2_Multiply_CCA_Axes,"_cca");
				testprop=makeName(options.Prop2_Multiply_CCA_Axes,"");
				testData1=getfield(csp1.xval_Dataset{i},testprop{:});
				csp1.xval_Dataset{i}=setfield(csp1.xval_Dataset{i},testccaprop{:},pagemtimes(testData1,W1));
				testData2=getfield(csp2.xval_Dataset{i},testprop{:});
				csp2.xval_Dataset{i}=setfield(csp2.xval_Dataset{i},testccaprop{:},pagemtimes(testData2,W2));
				% calculate correlation
				testccaCorrprop=makeName(options.Prop2_Multiply_CCA_Axes,"_corr");
				rho=diag(corr(f(testData1),f(testData2)));
				csp1.xval_Dataset{i}=setfield(csp1.xval_Dataset{i},testccaCorrprop{:},rho);
				csp2.xval_Dataset{i}=setfield(csp2.xval_Dataset{i},testccaCorrprop{:},rho);
			end
		end
		
		function name=makeName(propStrArr,name2append)
			props=cellstr(propStrArr);
			props=cellfun(@string,props,UniformOutput=false);
			name=props;
			name{end}=name{end}+name2append;
		end
	end
end

%% GCCA
methods (Static)
	function gcca_xval_Dataset(csp_struct,options)
		arguments
			csp_struct	% struct, each field containing CSP
			options.Prop2ApplyCCA (1,:) string
			options.Prop2_Multiply_CCA_Axes (1,:) string=[]
			options.ObtainAnchor (1,1) logical =true	% obtain anchor for testing data as well
		end
		prop=makeName(options.Prop2ApplyCCA,"");
		ccaprop=makeName(options.Prop2ApplyCCA,"_cca");
		axprop=makeName(options.Prop2ApplyCCA,"_ccaAxes");
		corrprop=makeName(options.Prop2ApplyCCA,"_ccaCorr");
		anchorprop=makeName(options.Prop2ApplyCCA,"_ccaAnchor");

		flds=string(fieldnames(csp_struct));
		csp_cell=struct2cell(csp_struct);		
		for i=1:numel(csp_struct.(flds(1)).xval_Dataset)	% all cross validation groups
			datacell=cellfun(@(csp)getfield(csp.xval_Dataset{i},prop{:}),csp_cell,UniformOutput=false);
			orgdatacell=datacell;
			f=@(a)multiBlockCat(a,[3:ndims(a)],1);
			datacell=cellfun(f,datacell,UniformOutput=false);
			GCCA=gcca(datacell{:});
			% reshape GCCA.Y to the usual block (assuming all the arrays
			% inside datacell share the same size)
			[flatdata,info]=multiBlockCat(orgdatacell{1},[3:ndims(orgdatacell{1})],1);
			if ~isequal(size(flatdata),size(GCCA.Y))
				error("size of flatdata and GCCA.Y must match");
			end
			GCCA.Y=multiBlockCatRevert(GCCA.Y,info);	% now GCCA.Y is in block form, corresponding to time x cell x trial x direction

			for s_ind=1:numel(flds)
				s=flds(s_ind);
				csp=csp_struct.(s);
				csp.xval_Dataset{i}=setfield(csp.xval_Dataset{i},axprop{:},GCCA.A{s_ind});
				csp.xval_Dataset{i}=setfield(csp.xval_Dataset{i},ccaprop{:},pagemtimes(orgdatacell{s_ind}-GCCA.mean{s_ind},GCCA.A{s_ind}));	% subtract mean for proper projection compatible with anchor
				csp.xval_Dataset{i}=setfield(csp.xval_Dataset{i},corrprop{:},GCCA.R{s_ind});
				csp.xval_Dataset{i}=setfield(csp.xval_Dataset{i},anchorprop{:},GCCA.Y);
				if ~isempty(options.Prop2_Multiply_CCA_Axes)
					testccaprop=makeName(options.Prop2_Multiply_CCA_Axes,"_cca");
					testprop=makeName(options.Prop2_Multiply_CCA_Axes,"");
					testData=getfield(csp.xval_Dataset{i},testprop{:});
					csp.xval_Dataset{i}=setfield(csp.xval_Dataset{i},testccaprop{:},pagemtimes(testData-GCCA.mean{s_ind},GCCA.A{s_ind}));	% subtract mean for proper projection compatible with anchor
				end
			end
			if options.ObtainAnchor
				XW_i={};
				for s_ind=1:numel(flds)
					s=flds(s_ind);
					XiWi=getfield(csp_struct.(s).xval_Dataset{i},testccaprop{:});	% (X_i-mean(trainX_i,1))*W_i where X_i is the test data and W_i is the weight for gcca, keep in mind that X_i is a multidimensional array
					XW_i{s_ind}=XiWi;
				end
				maxdim=max(cellfun(@(A)ndims(A),XW_i,UniformOutput=true))+1;
				XW_cat=cat(maxdim,XW_i{:});
				XW_mean=mean(XW_cat,maxdim);	% theoretically, the anchor y=sum(X_i*w_i, in terms of i) --> Y=sum(X_i*W_i) when Y and W_i have columns corresponding to cca components
				testanchorprop=makeName(options.Prop2_Multiply_CCA_Axes,"_ccaAnchor");
				testccaCorrprop=makeName(options.Prop2_Multiply_CCA_Axes,"_ccaCorr");
				for s_ind=1:numel(flds)
					s=flds(s_ind);
					csp=csp_struct.(s);
					csp.xval_Dataset{i}=setfield(csp.xval_Dataset{i},testanchorprop{:},XW_mean);
					XiWi=getfield(csp_struct.(s).xval_Dataset{i},testccaprop{:});
					rho=diag(corr(f(XiWi),f(XW_mean)))';	% correlation between test data canonical components and anchors
					csp.xval_Dataset{i}=setfield(csp.xval_Dataset{i},testccaCorrprop{:},rho);
				end
			end
		end
		
		function name=makeName(propStrArr,name2append)
			props=cellstr(propStrArr);
			props=cellfun(@string,props,UniformOutput=false);
			name=props;
			name{end}=name{end}+name2append;
		end
	end
end


%% TreeBagger decoding %%%%%%%%%%%%%%%%%%%%%%%%%%
methods
	function TreeBaggerDecode(obj,options)
		arguments
			obj
			options.NumTrees=100;
			options.Predictor2Use (1,1) string
			options.Responser2Use (1,1) string
			options.UseParallel (1,1) logical =false
			options.Regression (1,1) logical =false
			options.DisplayTime (1,1) logical =true
			options.DisplayText (1,1) logical =true
		end
		regressOn=options.Regression;
		NumTrees=options.NumTrees;
		if options.UseParallel
			paroptions = statset('UseParallel', true);
		end
		
		for i=1:numel(obj.xval_Dataset)
			tic
			Trainer=obj.xval_Dataset{i}.Trainer;
			Tester=obj.xval_Dataset{i}.Tester;
			pred=Trainer.(options.Predictor2Use);
			resp=Trainer.(options.Responser2Use);
			f=@(a)multiBlockCat(a,[3:ndims(a)],1);
			pred=f(pred);
			resp=f(resp);
			if size(pred, 1) ~= size(resp, 1)
				error('The number of timepoints in pred and resp must match.');
			end
			test_pred=Tester.(options.Predictor2Use);
			test_resp=Tester.(options.Responser2Use);
			flat_test_pred=f(test_pred);
			[flat_test_resp,info]=multiBlockCat(test_resp,[3:ndims(test_resp)],1);
			flat_rec_resp=zeros(size(flat_test_resp));
			if options.UseParallel
				for j=1:size(resp,2)
					if regressOn
						y=resp(:,j);
						Mdl=TreeBagger(NumTrees, ...
						pred,y, ...
						Method="regression", ...
						OOBPrediction="off",...
						Options=paroptions);
						flat_rec_resp(:,j)=Mdl.predict(flat_test_pred);
					else
						y=categorical(resp(:,j));
						Mdl=TreeBagger(NumTrees, ...
						pred,y, ...
						Method="classification", ...
						OOBPrediction="off",...
						Options=paroptions);
						flat_rec_resp(:,j)=str2double(Mdl.predict(flat_test_pred));
					end
				end
			else
				parfor j=1:size(resp,2)
					if regressOn
						y=resp(:,j);
						Mdl=TreeBagger(NumTrees, ...
						pred,y, ...
						Method="regression", ...
						OOBPrediction="off");
						flat_rec_resp(:,j)=Mdl.predict(flat_test_pred);
					else
						y=categorical(resp(:,j));
						Mdl=TreeBagger(NumTrees, ...
						pred,y, ...
						Method="classification", ...
						OOBPrediction="off");
						flat_rec_resp(:,j)=str2double(Mdl.predict(flat_test_pred));
					end
				end
			end
			if ~isequal(size(flat_rec_resp),size(flat_test_resp))
				error("The size of original and reconstructed data must match")
			end
			rec_resp=multiBlockCatRevert(flat_rec_resp,info);
			if ~isequal(size(rec_resp),size(test_resp))
				error("The size of original and reconstructed data must match")
			end
			obj.xval_Dataset{i}.Tester.Reconstructed=rec_resp;
			if options.DisplayText
				fprintf("xvalidation group %d / %d\n",i,numel(obj.xval_Dataset));
			end
			if options.DisplayTime
				toc
			end
		end
	end
end
%% ridge regression
methods
	function Weights=RidgeDecode(obj,options)
		arguments
			obj
			options.ridgeparams
			options.make_Data_single (1,1) logical =false
			options.UseGPU (1,1) logical =true
			options.Predictor2Use (1,1) string	% Trainer.(options.Predictor2Use)
			options.Responser2Use (1,1) string	% Trainer.(options.Responser2Use);
		end
		obj.ridgeparams=options.ridgeparams;
		meanWeights=[];
		for i=1:numel(obj.xval_Dataset)
			tic
			Trainer=obj.xval_Dataset{i}.Trainer;
			Tester=obj.xval_Dataset{i}.Tester;
			pred=Trainer.(options.Predictor2Use);
			resp=Trainer.(options.Responser2Use);
			f=@(a)multiBlockCat(a,[3:ndims(a)],1);
			pred=f(pred);
			resp=f(resp);
			if size(pred, 1) ~= size(resp, 1)
				error('The number of timepoints in pred and resp must match.');
			end
			test_pred=Tester.(options.Predictor2Use);
			test_resp=Tester.(options.Responser2Use);
			flat_test_pred=f(test_pred);
			[flat_test_resp,info]=multiBlockCat(test_resp,[3:ndims(test_resp)],1);
			% training
			if options.make_Data_single
				x = (single(pred));
				y = (single(resp));
			else
				x = (pred);
				y = (resp);
			end
			if options.UseGPU
				x=gpuArray(x);
				y=gpuArray(y);
			end
			B = ridge_optimized(y, x, options.ridgeparams);	
			clear y x;

			% Process tester data
			if options.make_Data_single
				x = (single(flat_test_pred));
				y = (single(flat_test_resp));
			else
				x = (flat_test_pred);
				y = (flat_test_resp);
			end
			if options.UseGPU
				x=gpuArray(x);
				y=gpuArray(y);
			end
			y_recon = B(1,:,:) + pagemtimes(x, B(2:end,:,:));

			if options.UseGPU
				B=gather(B);
				y=gather(y);
				y_recon=gather(y_recon);	% time x parameters x ridgeparams
			end
			if options.make_Data_single
				B=double(B);
				y=double(y);
				y_recon=double(y_recon);
			end
			rec_ridgeparams={};
			for ridgepara=1:size(y_recon,3)
				flat_rec_resp=y_recon(:,:,ridgepara);
				if ~isequal(size(flat_rec_resp),size(flat_test_resp))
					error("The size of original and reconstructed data must match")
				end
				rec_resp=multiBlockCatRevert(flat_rec_resp,info);
				if ~isequal(size(rec_resp),size(test_resp))
					error("The size of original and reconstructed data must match")
				end
				rec_ridgeparams{ridgepara}=rec_resp;
			end
			obj.xval_Dataset{i}.Tester.Reconstructed=rec_ridgeparams;	% cell array with element --> time x params x trials...

			% get weights
			if i==1
				meanWeights=B;
			else
				meanWeights=(1/i)*(meanWeights.*(i-1)+B);
			end

			fprintf("xvalidation group %d / %d\n",i,numel(obj.xval_Dataset));
			toc
		end
		% get weights %%%%%%%%%%%
		Weights=meanWeights;
	end
	function R2_ridge_cells = xval_assess_R2_ridge(obj, options)
		% R2_ridge_cells: ridgeparams x response-parameters
		%
		% This reproduces the same R used in xval_assess_performance():
		%   R = 1 - sum((rec-resp).^2,1) ./ sum((resp-resp_null).^2,1)
		% where resp_null is the mean of the TRAINING fold, evaluated on the TEST fold.
		%
		% Unlike xval_assess_performance_ridge(), this method:
		%   - does not copy the object
		%   - does not concatenate all folds
		%   - computes only R^2

		arguments
			obj
			options.ReferenceResponser (1,1) string
		end

		nRidge = numel(obj.ridgeparams);

		SSE = [];   % ridgeparams x response-parameters
		SST0 = [];  % 1 x response-parameters

		for i = 1:numel(obj.xval_Dataset)
			train_resp = obj.xval_Dataset{i}.Trainer.(options.ReferenceResponser);
			test_resp  = obj.xval_Dataset{i}.Tester.(options.ReferenceResponser);
			rec_ridgeparams = obj.xval_Dataset{i}.Tester.Reconstructed;

			if ~iscell(rec_ridgeparams)
				error("Tester.Reconstructed must be a cell array with one element per ridge parameter.")
			end
			if numel(rec_ridgeparams) ~= nRidge
				error("numel(Tester.Reconstructed) must match numel(obj.ridgeparams).")
			end

			% Mean of the training fold, preserving parameter dimension (dim 2)
			mu_train = meanOverObservationDims(train_resp);              % 1 x params x 1 x ...

			% Denominator for this fold: sum over all observations in the TEST fold
			sst0_i = double(sumSqOverObservationDims(test_resp - mu_train));   % 1 x params

			if isempty(SSE)
				nResp = numel(sst0_i);
				SSE  = zeros(nRidge, nResp);
				SST0 = zeros(1, nResp);
			else
				if numel(sst0_i) ~= size(SSE, 2)
					error("The number of response parameters changed across cross-validation folds.")
				end
			end

			SST0 = SST0 + sst0_i;

			for r = 1:nRidge
				rec_resp = rec_ridgeparams{r};

				if ~isequal(size(rec_resp), size(test_resp))
					error("The size of original and reconstructed data must match.")
				end

				SSE(r, :) = SSE(r, :) + double(sumSqOverObservationDims(rec_resp - test_resp));
			end

			fprintf("R2 assessment group %d / %d\n", i, numel(obj.xval_Dataset));
		end

		R2_ridge_cells = 1 - (SSE ./ SST0);

		function mu = meanOverObservationDims(a)
			% Average over all observation dimensions, keeping dim 2 (= parameters)
			mu = mean(a, 1);
			for d = 3:ndims(a)
				mu = mean(mu, d);
			end
		end

		function s = sumSqOverObservationDims(a)
			% Sum squared values over all observation dimensions, keeping dim 2
			s = sum(a.^2, 1);
			for d = 3:ndims(a)
				s = sum(s, d);
			end
			s = reshape(s, 1, []);
		end
	end
	function Scores_ridge_cell=xval_assess_performance_ridge(obj,options)
		% Score_ridge_cell: struct with each field being ridge params x
		% neurons
		arguments
			obj
			options.ReferenceResponser (1,1) string
		end
		csp=copy(obj);
		scores_ridgeparams={};
		for r=1:numel(obj.ridgeparams)
			Csp=copy(obj);
			for i=1:numel(obj.xval_Dataset)
				Csp.xval_Dataset{i}.Tester.Reconstructed=obj.xval_Dataset{i}.Tester.Reconstructed{r};
			end
			scores_ridgeparams{r}=Csp.xval_assess_performance(ReferenceResponser=options.ReferenceResponser);
		end
		score_type_ridgeparams=CellWithStruct2StructWithCell(scores_ridgeparams);
		Scores_ridge_cell=organizeScores(score_type_ridgeparams);
		function scores_ridge_cells=organizeScores(scores)
			scores_ridge_cells=[];
			flds=string(fields(scores));
			for s=flds(:)'
				scores_ridge_cells.(s)=cat(1,scores.(s){:});
			end
		end
	end
	function R_null = xval_assess_testmean_null_R(obj, options)
		arguments
			obj
			options.ReferenceResponser (1,1) string
		end

		resp_g = {};
		rec_null_g = {};
		resp_null_train_g = {};

		f = @(a) multiBlockCat(a, [3:ndims(a)], 1);

		for i = 1:numel(obj.xval_Dataset)
			train_resp = f(obj.xval_Dataset{i}.Trainer.(options.ReferenceResponser));
			test_resp  = f(obj.xval_Dataset{i}.Tester.(options.ReferenceResponser));

			train_mean = mean(train_resp, 1);
			test_mean  = mean(test_resp, 1);

			resp_g{i} = test_resp;
			rec_null_g{i} = zeros(size(test_resp), 'like', test_resp) + test_mean;
			resp_null_train_g{i} = zeros(size(test_resp), 'like', test_resp) + train_mean;
		end

		resp = cat(1, resp_g{:});
		rec_null = cat(1, rec_null_g{:});
		resp_null_train = cat(1, resp_null_train_g{:});

		denom = sum((resp - resp_null_train).^2, 1);
		numer = sum((rec_null - resp).^2, 1);

		R_null = nan(size(denom));
		valid = denom > 0;
		R_null(valid) = 1 - (numer(valid) ./ denom(valid));
	end
end
%%
properties

end
methods
	function Scores=xval_assess_performance(obj,options)	% Tester.Reconstructed vs Tester.(options.ReferenceResponser)
		arguments
			obj
			options.ReferenceResponser (1,1) string
		end
		% store all data in cell arrays
		train_resp_g={};
		resp_g={};
		rec_g={};
		f=@(a)multiBlockCat(a,[3:ndims(a)],1);
		for i=1:numel(obj.xval_Dataset)
			tic
			train_resp_g{i}=f(obj.xval_Dataset{i}.Trainer.(options.ReferenceResponser));
			resp_g{i}=f(obj.xval_Dataset{i}.Tester.(options.ReferenceResponser));
			rec_g{i}=f(obj.xval_Dataset{i}.Tester.Reconstructed);
		end
		test_null_g=make_null_model(train_resp_g,resp_g);

		% concatenate all data
		resp=cat(1,resp_g{:});
		rec=cat(1,rec_g{:});
		resp_null=cat(1,test_null_g{:});

		Scores=[];
		Scores.expVar=ExplainedVariance(resp_g,rec_g);
		Scores.R=Rsquared(resp,rec,resp_null);
		Scores.MeanSquaredError=MeanSquaredError(resp,rec);
		Scores.MeanAbsoluteError=MeanAbsoluteError(resp,rec);
		Scores.PearsonCorrelation=PearsonCorrelation(resp_g,rec_g);

		function test_null_g=make_null_model(train_resp_g,resp_g)
			test_null_g=cellfun(@(train_resp,test_resp)(zeros(size(test_resp))+mean(train_resp,1)), ...
				train_resp_g,resp_g,UniformOutput=false);
		end
		function expVar=ExplainedVariance(resp_g,rec_g)
			score=[];
			for g=1:numel(resp_g)
				y=resp_g{g};
				y_recon=rec_g{g};
				score(g,:) = 1 - (var(y-y_recon,1,1)./var(y,1,1));
			end
			expVar=mean(score,1);
		end
		function R=Rsquared(resp,rec,resp_null)
			R=1-(sum((rec-resp).^2,1)./sum((resp-resp_null).^2));
		end
		function score = MeanSquaredError(resp,rec)
			y=resp;
			y_recon=rec;
			score = sum((y - y_recon).^2, 1) / size(y, 1);
		end
		function score = MeanAbsoluteError(resp,rec)
			y=resp;
			y_recon=rec;
			score = sum(abs(y - y_recon), 1) / size(y, 1);
		end
		function Corr = PearsonCorrelation(resp_g,rec_g)
			score=[];
			for g=1:numel(resp_g)
				y=resp_g{g};
				y_recon=rec_g{g};
				score(g,:) = fastPearsonCorr(y, y_recon,1);
			end
			Corr=mean(score,1);
		end
	end
end

%% analyze predicted xval_Dataset %%%%%%%%%%%%%%%%%%%%%%%
methods
	function Undiscretize_reconstructed_values(obj)
		for i=1:numel(obj.xval_Dataset)
			rec=obj.xval_Dataset{i}.Tester.Reconstructed;
			rec_cell=mat2cell_2ndDim(rec);
			funcs=obj.Discretized_SP_and_DLC_undiscretizer;
			if numel(funcs)~=numel(rec_cell)
				error("funcs and rec_cell must have the same size")
			end
			undis_cell=cellfun(@(a,b)b(a),rec_cell,funcs,UniformOutput=false);
			% undis_cell={};
			% for j=1:numel(funcs)
			% 	f=funcs{j};
			% 	undis_cell{j}=f(rec_cell{j});
			% end
			obj.xval_Dataset{i}.Tester.Rec_cell=rec_cell;
			obj.xval_Dataset{i}.Tester.Rec_cell_undiscretized=undis_cell;
			obj.xval_Dataset{i}.Tester.Resp_cell=mat2cell_2ndDim(obj.xval_Dataset{i}.Tester.Responser);
			obj.xval_Dataset{i}.Tester.Resp_cell_undiscretized=cellfun(@(a,b)b(a), ...
				obj.xval_Dataset{i}.Tester.Resp_cell,funcs,UniformOutput=false);
		end

		function newdata=mat2cell_2ndDim(data)
			currsize=num2cell(size(data));
			newdata=mat2cell(data,currsize{1},ones(1,currsize{2}),currsize{3:end});
		end
	end
end

%% ridge regression of xval_Dataset
properties	
	ridge_weights	% cell of groups
	ridge_y			% cell of groups
	ridge_y_rec		% cell of groups
	ridge_performance
	ridge_params
end
methods
	function RidgeRegression(obj,ridgeparams,options)
			arguments
				obj
				ridgeparams
				options.MakeDataSingle=false
				options.DisplayGPUMemory=false
				options.useGPU (1,1) =true

				options.Predictor2Use (1,1) string
				options.Responser2Use (1,1) string
			end
			Trainer=obj.Trainer;
			Tester=obj.Tester;
			B_group=[];
			y_group=[];
			y_recon_group=[];
			y_null_group=[];
			%%%%%%%%%%%%%%%%%%
			for i=1:numel(obj.xval_Dataset)
				Trainer=obj.xval_Dataset{i}.Trainer;
				Tester=obj.xval_Dataset{i}.Tester;

				train_pred=Trainer.(options.Predictor2Use);
				train_resp=Trainer.(options.Responser2Use);
				test_pred=Tester.(options.Predictor2Use);
				test_resp=Tester.(options.Responser2Use);
				if options.DisplayGPUMemory
					fprintf("Available GPU Memory: %d\n", gpuDevice().AvailableMemory);
				end
				 [B_group{i},y_group{i},y_recon_group{i},y_null_group{i}] = processGroup(train_pred,train_resp,test_pred,test_resp, ridgeparams, options.MakeDataSingle, i,options);
			end
			obj.ridge_params=ridgeparams;

			Meanscore=assessPerformance(y_group,y_recon_group,y_null_group);

			obj.ridge_performance=Meanscore;
			obj.ridge_weights=B_group;
			obj.ridge_y=y_group;
			obj.ridge_y_rec=y_recon_group;
			
			function [B,y,y_recon,y_null] = processGroup(train_pred,train_resp,test_pred,test_resp, ridgeparams, MakeDataSingle, groupIndex,options)
				% Process trainer data
				y=train_resp;
				x=train_pred;
				if MakeDataSingle
					y=single(y);
					x=single(x);
				end
				if options.useGPU
					y=gpuArray(y);
					x=gpuArray(x);
				end
				B = ridge_optimized(y, x, ridgeparams);
				y_train=gather(y);
				clear y x;

				% Process tester data
				y=test_resp;
				x=test_pred;
				if MakeDataSingle
					y = (single(y));
					x = (single(x));
				end
				if options.useGPU
					y=gpuArray(y);
					x=gpuArray(x);
				end

				y_recon = B(1,:,:) + pagemtimes(x, B(2:end,:,:));

				B=gather(B);
				y=gather(y);
				y_recon=gather(y_recon);				%% testing time points x cell x ridgeparams
				y=y+zeros(size(y_recon));				%% testing time points x cell x ridgeparams
				y_null=make_null_model(y_recon,y_train);%% testing time points x cell x ridgeparams
			end

			function y_train_mean=make_null_model(y_recon,y_train)
				y_train_mean=mean(y_train,1);					%% training time points x cell
				y_train_mean=y_train_mean+zeros(size(y_recon));	%% testing time points x cell x ridgeparams
			end

			function Meanscore=assessPerformance(y_group,y_recon_group,y_null_group)
				% arrays inside are cell x ridgeparams
				y=cell2mat(y_group(:));
				y_recon=cell2mat(y_recon_group(:));
				y_null=cell2mat(y_null_group(:));
				Meanscore.MeanSquaredError=MeanSquaredError(y, y_recon);
				Meanscore.MeanAbsoluteError=MeanAbsoluteError(y, y_recon);
				Meanscore.Rsquared=Rsquared(y, y_recon);
				Meanscore.pseudoRsquared=pseudoRsquared(y, y_recon,y_null);
				Meanscore.ExplainedVariance=ExplainedVariance(y, y_recon);
				Meanscore.PearsonCorrelation=PearsonCorrelation(y, y_recon);
			end

			
			% Below, first dim is time, second is cell, third is
			% ridgeparams
			% --> result: cell x ridgeparams
			function score = MeanSquaredError(y, y_recon)
				score = sum((y - y_recon).^2, 1) / size(y, 1);
				score=DeleteSingleDim(score,1);
			end
			function score = MeanAbsoluteError(y, y_recon)
				score = sum(abs(y - y_recon), 1) / size(y, 1);
				score=DeleteSingleDim(score,1);
			end
			function score = Rsquared(y, y_recon)
				normy = y - mean(y, 1);
				score = 1 - ( sum((y - y_recon).^2, 1) ./ sum(normy.^2, 1) );
				score=DeleteSingleDim(score,1);
			end
			function score = pseudoRsquared(y, y_recon,y_null)
				score=1-(sum((y-y_recon).^2,1)./sum((y-y_null).^2,1));	% 1 x cell x ridge params
				score=DeleteSingleDim(score,1);
			end
			function score = ExplainedVariance(y, y_recon)
				score=1-(var(y-y_recon,0,1)./var(y,0,1));
				score=DeleteSingleDim(score,1);
			end
			function score = PearsonCorrelation(y, y_recon)
				score = fastPearsonCorr(y, y_recon,1);
				score=DeleteSingleDim(score,1);
			end
		end
end

%% clear and organize xval_Dataset
properties
	xval_Dataset_cleared
end
methods
	function organize_xval_Dataset(obj,fields2keep)
		arguments
			obj
		end
		arguments (Repeating)
			fields2keep (1,:) string
		end
		obj.xval_Dataset_cleared=cell(size(obj.xval_Dataset));
		for i=1:numel(obj.xval_Dataset)
			for j=1:numel(fields2keep)
				currstr=fields2keep{j};
				currstr=convert2StrCell(currstr);
				data=getfield(obj.xval_Dataset{i},currstr{:});
				obj.xval_Dataset_cleared{i}=setfield(obj.xval_Dataset_cleared{i},currstr{:},data);
			end
		end

		function name=convert2StrCell(propStrArr)
			props=cellstr(propStrArr);
			name=cellfun(@string,props,UniformOutput=false);
		end
	end
	%% Legacy (same as organize_xval_Dataset but the name is confusing)
	function clear_xval_Dataset(obj,fields2keep)
		arguments
			obj
		end
		arguments (Repeating)
			fields2keep (1,:) string
		end
		obj.xval_Dataset_cleared=cell(size(obj.xval_Dataset));
		for i=1:numel(obj.xval_Dataset)
			for j=1:numel(fields2keep)
				currstr=fields2keep{j};
				currstr=convert2StrCell(currstr);
				data=getfield(obj.xval_Dataset{i},currstr{:});
				obj.xval_Dataset_cleared{i}=setfield(obj.xval_Dataset_cleared{i},currstr{:},data);
			end
		end

		function name=convert2StrCell(propStrArr)
			props=cellstr(propStrArr);
			name=cellfun(@string,props,UniformOutput=false);
		end
	end
end



%% split CSP to xvalidation groups
properties
	cspArray
	DataSet
end
methods
	function cspArray=split_CSP_to_xvalidation_groups(obj,options)
		arguments
			obj
			options.groupnum
		end
		% Trainer=cell(1,options.groupnum);
		% Tester=cell(1,options.groupnum);
		Trainer=[];
		Tester=[];
		dim2split=3;
		[Trainer,Tester]=splitArray(Trainer,Tester,"tck",obj.tck,dim2split);
		[Trainer,Tester]=splitArray(Trainer,Tester,"tpk",obj.tpk,dim2split);

		% split DLC data
		Trainer.DLC=struct;
		Tester.DLC=struct;
		numDLC=obj.DLC_convert2numeric(obj.DLC);
		flds=string(fields(numDLC));
		for fld=flds(:)'
			[Trainer.DLC.(fld),Tester.DLC.(fld)]=splitArrayToGroups(numDLC.(fld),options.groupnum,dim2split);
			Trainer.DLC=csp.DLC_convert2cell(Trainer.DLC);
			Tester.DLC=csp.DLC_convert2cell(Tester.DLC);
		end

		DataSet=[];
		DataSet.Trainer=Trainer;
		DataSet.Tester=Tester;

		DataSet=StructWithCell2CellWithStruct_multiLayers(DataSet);

		TrainerCSP=CSP().empty(0,length(DataSet));
		TesterCSP=CSP().empty(0,length(DataSet));
		for g=1:options.groupnum
			TrainerCSP(g)=copy(csp);
			TesterCSP(g)=copy(csp);
			Trainer=DataSet{g}.Trainer;
			Tester=DataSet{g}.Tester;
			flds=string(fields(Trainer));
			for fld=flds(:)'
				TrainerCSP(g).(fld)=Trainer.(fld);
			end
			flds=string(fields(Tester));
			for fld=flds(:)'
				TesterCSP(g).(fld)=Tester.(fld);
			end
		end
		cspArray=[];
		cspArray.Trainer=TrainerCSP;
		cspArray.Tester=TesterCSP;
		obj.cspArray=cspArray;

		function [Trainer,Tester]=splitArray(Trainer,Tester,fieldName,array,dim2split)
			[Trainer.(fieldName),Tester.(fieldName)]=splitArrayToGroups(array, options.groupnum, dim2split);
		end
	end
	function split_cell_sp_dlc_to_xvalidation_groups(obj,options)
		arguments
			obj
			options.groupnum
		end
		f=@(a)multiBlockCat(a,[3:ndims(a)],1);
		numDLC=obj.DLC_convert2numeric(obj.DLC);
		tckd=obj.tck;
		tpkd=obj.tpk;
		Trainer=[];
		Tester=[];
		dim2split=3;
		[Trainer.TKDc,Tester.TKDc]=splitArrayToGroups(tckd, options.groupnum, dim2split);
		[Trainer.TKDp,Tester.TKDp]=splitArrayToGroups(tpkd, options.groupnum, dim2split);
		Trainer.DLC=struct;
		Tester.DLC=struct;
		flds=string(fields(numDLC));
		for fld=flds(:)'
			[Trainer.DLC.(fld),Tester.DLC.(fld)]=splitArrayToGroups(numDLC.(fld),options.groupnum,dim2split);
		end
		obj.DataSet.Trainer=Trainer;
		obj.DataSet.Tester=Tester;
		obj.DataSet=applyFunctionLayers(obj.DataSet,f);
		function convertStructWithCell_to_CellWithStruct(cellArray,Struct,fieldname)
			Struct.(fieldname)={};
			for i=1:numel(cellArray)
				Struct.(fieldname){i}=cellArray{i};
			end
		end
	end
end
% %% dimensionality reduction
% properties
% 	tckd_pca
% 	tckd_pca_axes
% end
% methods
% 	function pca_tckd(obj,options)
% 		arguments
% 			obj
% 			options.NumDim2keep=[]
% 		end
% 		f=@(a)multiBlockCat(a,[3:ndims(a)],1);
% 		TKDc=f(obj.tck);
% 		Axes=pca(TKDc);
% 		numdim=size(Axes,2);
% 		if ~isempty(options.NumDim2keep)
% 			numdim=options.NumDim2keep;
% 		end
% 		IND=makeIND(Axes,2,1:numdim);
% 		obj.tckd_pca_axes=Axes(IND{:});
% 		obj.tckd_pca=pagemtimes(obj.tck,Axes(IND{:}));
% 	end
% 	% 
% 	% function pca_on_DataSet_TKDc(obj,options)
% 	% 	arguments
% 		% 	obj
% 		% 	options.NumDim2keep=[]
% 	% 	end
% 	% 	obj.DataSet.Trainer.TKDc_pca_axes={};
% 	% 	obj.DataSet.Trainer.TKDc_pca={};
% 	% 	obj.DataSet.Tester.TKDc_pca={};
% 	% 	for i=1:numel(obj.DataSet.Trainer.TKDc)
% 		% 	Axes=pca(obj.DataSet.Trainer.TKDc{i});
% 		% 	numdim=size(Axes,2);
% 		% 	if ~isempty(options.NumDim2keep)
% 			% 	numdim=options.NumDim2keep;
% 		% 	end
% 		% 	IND=makeIND(Axes,2,1:numdim);
% 		% 	obj.DataSet.Trainer.TKDc_pca_axes{i}=Axes(IND{:});
% 		% 	obj.DataSet.Trainer.TKDc_pca{i}=obj.DataSet.Trainer.TKDc{i}*Axes(IND{:});
% 		% 	obj.DataSet.Tester.TKDc_pca{i}=obj.DataSet.Tester.TKDc{i}*Axes(IND{:});
% 	% 	end
% 	% end
% end
% properties
% 	tckd_cca
% end
% methods
% 	function cca_tckd(obj,csp1,csp2,cell_field2use)
% 
% 
% 	end
% end

%% ridge regression to decode posture and sp
properties
	combined_posture_4_decode
	posture_4_decode_params
end
methods
	function make_posture_for_decode(obj,params2decode,options)
		arguments
			obj
		end
		arguments (Repeating)
			params2decode (1,1) string {mustBeMember(params2decode,["tXYZbp","tDXYZbp","tD2XYZbp", ...
					"tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp", ...
					"tDXYZbp_magnitude","tDDOTbp_magnitude","tDANGbp_magnitude", ...
					"tD2XYZbp_magnitude","tD2DOTbp_magnitude","tD2ANGbp_magnitude"])}
		end
		arguments
			options.take_zscore (1,1) logical =true
		end
		params2decode=string(params2decode);
		numDLC=obj.DLC_convert2numeric(obj.DLC);
		XYZ=["tXYZbp","tDXYZbp","tD2XYZbp"];
		ANG=["tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"];
		speed=["tDXYZbp_magnitude","tDDOTbp_magnitude","tDANGbp_magnitude"];
		accel=["tD2XYZbp_magnitude","tD2DOTbp_magnitude","tD2ANGbp_magnitude"];
		
		combined_posture={};
		ind_2ndDim_for_each_posture={};
		accum_ind=1;
		for i=1:length(params2decode)
			param=params2decode(i);
			if any(XYZ==param)||any(ANG==param)
				data=numDLC.(param);
			elseif any(speed==param)
				data=obj.DLC_speed.(param);
			elseif any(accel==param)
				data=obj.DLC_accel.(param);
			else
				error("No specified data found, check params2decode")
			end
			if options.take_zscore
				data=zscore(data,0,[1,3:ndims(data)]);
			end
			combined_posture{end+1}=data;
			ind_2ndDim_for_each_posture{end+1}=accum_ind:(accum_ind+size(data,2)-1);
			accum_ind=accum_ind+size(data,2);
		end
		obj.combined_posture_4_decode=combined_posture;
		obj.posture_4_decode_params=[];
		obj.posture_4_decode_params.posture_ind=params2decode;
		obj.posture_4_decode_params.ind_2ndDim_for_each_posture=ind_2ndDim_for_each_posture;
	end
end
properties
	tck_4_decode
	cell_reg_4_decode
end
methods
	%% preprocess cell data for decoding
	function make_cellac_by_regions_for_decode(obj,options)
		arguments
			obj
			options.region_4_decode
			options.cell_num_2_randomly_sample=[]
			options.skew_thresh=[]
		end
		obj.cell_reg_4_decode=options.region_4_decode;
		reg_ind=obj.Reg_RegionNames==options.region_4_decode;
		reg_log=obj.Reg_inRegLog_reg{reg_ind};
		if isempty(options.skew_thresh)
			all_log=reg_log;
		else
			skew_log=obj.Skew>=options.skew_thresh;
			all_log=logical(reg_log(:).*skew_log(:));
		end
		if ~isempty(options.cell_num_2_randomly_sample)
			if sum(all_log)<options.cell_num_2_randomly_sample
				warning(sprintf("%s contains few cells (%d), returning empty array", ...
					obj.cell_reg_4_decode,sum(all_log)));
				all_log=[];
			else
				all_log=find(all_log);
				ind=randperm(length(all_log),options.cell_num_2_randomly_sample);
				all_log=all_log(ind);
			end
		end

		ind=makeIND(obj.tck,2,all_log);
		tckd=obj.tck(ind{:});

		obj.tck_4_decode=tckd;
	end
	%% make cross validation group
	function make_xvalidation_group_for_decode(obj,options)
		arguments
			obj
			options.groupnum
			options.Pred=[]
			options.Resp=[]
		end
		obj.groupnum=options.groupnum;
		if isempty(options.Pred)
			Pred=obj.tck_4_decode;
		else
			Pred=options.Pred;
		end
		if isempty(options.Resp)
			Resp=cell2mat(obj.combined_posture_4_decode);
		else
			Resp=options.Resp;
		end
		[obj.Trainer,obj.Tester]=splitPredResp(Pred,Resp,options.groupnum,3);
		% concatenate all the 3rd and higher dimensions to 1st for decoding
		f=@(a)multiBlockCat(a,[3:ndims(a)],1);
		obj.Trainer=applyFunctionLayers(obj.Trainer,f);
		obj.Tester=applyFunctionLayers(obj.Tester,f);
	end
	%%
	function combined_posture_cell=revert_combined_posture_to_cell(obj,combined_posture_numeric,dim2split)
		combined_posture_cell={};
		for i=1:length(obj.posture_4_decode_params.ind_2ndDim_for_each_posture)
			ind=createIND(combined_posture_numeric,dim2split,obj.posture_4_decode_params.ind_2ndDim_for_each_posture{i});
			combined_posture_cell{i}=combined_posture_numeric(ind{:});
		end
		if false
			%% check if this function works correctly
			isequal(csp.combined_posture_4_decode,2, ...
				csp.revert_combined_posture_to_cell(cell2mat(csp.combined_posture_4_decode)));
		end
	end
end
%% train ridge decoder
properties
	MeanScore_RespRidgeparams
	StdScore_RespRidgeparams
end
methods
	function train_ridge_decoder(obj,options)
		arguments
			obj
			options.ridgeparams
			options.nameOfCurrentSetting (1,1) string
			options.make_Data_single (1,1) logical =false
		end
		obj.ridgeparams=options.ridgeparams;
		B_group={};
		y_group={};
		y_recon_group={};
		for i=1:obj.groupnum
			[B_group{i},y_group{i},y_recon_group{i}] = processGroup(obj.Trainer, obj.Tester, options.ridgeparams, options.make_Data_single, i);
		end
		[Meanscore,Stdscore]=calculateRegressionPerformance(y_group,y_recon_group);
		obj.MeanScore_RespRidgeparams=Meanscore;
		obj.StdScore_RespRidgeparams=Stdscore;
		
		function [B,y,y_recon] = processGroup(Trainer, Tester, ridgeparams, MakeDataSingle, groupIndex)
			% Process trainer data
			if MakeDataSingle
				y = gpuArray(single(Trainer.Resp{groupIndex}));
				x = gpuArray(single(Trainer.Pred{groupIndex}));
			else
				y = gpuArray(Trainer.Resp{groupIndex});
				x = gpuArray(Trainer.Pred{groupIndex});
			end
			B = ridge_optimized(y, x, ridgeparams);
			clear y x;

			% Process tester data
			if MakeDataSingle
				x = gpuArray(single(Tester.Pred{groupIndex}));
				y = gpuArray(single(Tester.Resp{groupIndex}));
			else
				x = gpuArray(Tester.Pred{groupIndex});
				y = gpuArray(Tester.Resp{groupIndex});
			end
			y_recon = B(1,:,:) + pagemtimes(x, B(2:end,:,:));

			B=gather(B);
			y=gather(y);
			y_recon=gather(y_recon);
		end
		function [Meanscore,Stdscore]=calculateRegressionPerformance(y_group,y_recon_group)
			y_group=permute(y_group(:),[2,3,4,1]);
			y_recon_group=permute(y_recon_group(:),[2,3,4,1]);
			Meanscore=[];
			Stdscore=[];

			[Meanscore.MeanSquaredError,Stdscore.MeanSquaredError]=AssessPerformanceAllCell(y_group,y_recon_group,@MeanSquaredError);
			[Meanscore.MeanAbsoluteError,Stdscore.MeanAbsoluteError]=AssessPerformanceAllCell(y_group,y_recon_group,@MeanAbsoluteError);
			[Meanscore.Rsquared,Stdscore.Rsquared]=AssessPerformanceAllCell(y_group,y_recon_group,@Rsquared);
			[Meanscore.ExplainedVariance,Stdscore.ExplainedVariance]=AssessPerformanceAllCell(y_group,y_recon_group,@ExplainedVariance);
			[Meanscore.PearsonCorrelation,Stdscore.PearsonCorrelation]=AssessPerformanceAllCell(y_group,y_recon_group,@PearsonCorrelation);
		end
		function [meanscore,stdscore]=AssessPerformanceAllCell(y_group,y_recon_group,func_handle)
			score_group=cellfun(func_handle,y_group,y_recon_group,UniformOutput=false);
			score_group=cell2mat(score_group);	% 1 x responser x ridgeparams x groups
			score_group=DeleteSingleDim(score_group,1);	% responser x ridgeparams x groups
			meanscore=mean(score_group,3); % responser x ridgeparams
			stdscore=std(score_group,1,3); % responser x ridgeparams
		end
		% Below, first dim is time, second is cell, third is
		% ridgeparams
		% --> result: 1 x cell x ridgeparams
		function score = MeanSquaredError(y, y_recon)
			score = sum((y - y_recon).^2, 1) / size(y, 1);
		end
		function score = MeanAbsoluteError(y, y_recon)
			score = sum(abs(y - y_recon), 1) / size(y, 1);
		end
		function score = Rsquared(y, y_recon)
			y = y + zeros(size(y_recon), class(y_recon));  % broadcast to match dimensions
			normy = y - mean(y, 1);
			score = 1 - ( sum((y - y_recon).^2, 1) ./ sum(normy.^2, 1) );
		end
		function score = ExplainedVariance(y, y_recon)
			y = y + zeros(size(y_recon), class(y_recon));  % broadcasting
			score = 1 - ( sum((y - y_recon).^2, 1) ./ sum(y.^2, 1) );
		end
		function score = PearsonCorrelation(y, y_recon)
			score = fastPearsonCorr(y, y_recon,1);
		end
	end
	
end

%% Eigen Posture %%%%%%%%%%%%%%%%%
properties
	eigen_DLC
	eigen_DLC_axes
	eigen_DLC_edges
	eigen_DLC_var
	eigen_DLC_mean
end
methods
	function make_eigen_posture(obj,options)
		arguments
			obj
			options.NumDim2keep=[]
			options.zscore_data (1,1) logical=false
		end
		numDLC=obj.DLC_convert2numeric(obj.DLC);
		if options.zscore_data
			numDLC.tXYZbp=zscore(numDLC.tXYZbp,0,[1,3:ndims(numDLC.tXYZbp)]);
		end
		f=@(a)multiBlockCat(a,[3:ndims(a)],1);
		flatDLC=f(numDLC.tXYZbp);
		[Axes,~,posvar]=pca(flatDLC);
		if ~isempty(options.NumDim2keep)
			Axes=Axes(:,1:options.NumDim2keep);
		end
		obj.eigen_DLC=pagemtimes(numDLC.tXYZbp,Axes);
		obj.eigen_DLC_axes=Axes;
		obj.eigen_DLC_var=posvar;
		obj.eigen_DLC_mean=mean(flatDLC,1);
	end
	function recX=reconstruct_from_eigen_posture(obj,options)
		arguments
			obj
			options.Dims2use (1,:)
			options.PostureData=[]
		end
		if isempty(options.PostureData)
			numDLC=obj.cat3Dcoordinates(obj.DLC.tXYZbp);
		else
			numDLC=obj.cat3Dcoordinates(options.PostureData);
		end
		f=@(a)multiBlockCat(a,[3:ndims(a)],1);
		X=f(numDLC);
		meanX=mean(X,1);
		f=@(a,b)pagemtimes(a,b);
		W=[obj.eigen_DLC_axes(:,options.Dims2use)];
		recX=f(f(numDLC-meanX,W),W')+meanX;
	end
end


%% decoding and encoding analysis
	properties
		SP_edges_cellVector
		SP_bincenter_cellVector
		SP_Discretized_cellVector
		SP_Discretized_names
		SP_Discretized_Discretizer
		SP_Discretized_Undiscretizer
		SP_Discretized_Undiscretized
		SP_Discretized_edges
		SP_Discretized_bincount
		SP_Discretized_params2encode
		SP_Discretized_params2encode_cell
		SP_Discretized_orgData
		SP_Discretized_Undiscretized_names

		DLC_inv_or_normal
		DLC_edges_type_cellVector
		DLC_bincenter_type_cellVector
		DLC_edges_speed_cellVector
		DLC_speed
		DLC_bincenter_speed_cellVector
		DLC_edges_speed_cellVector_names
		DLC_edges_accel_cellVector
		DLC_accel
		DLC_bincenter_accel_cellVector
		DLC_edges_accel_cellVector_names
		DLC_Discretized_bincount
		DLC_Discretized_edge_percent_thresh_For_cutoff

		Trainer
		Tester
		Rec		% Corresponds with Tester
		groupnum
	end
	methods
		function [tXYZ_speed,tANG_speed,tDOT_speed]=make_DLC_speed(obj)
			tXYZ_speed=cellfun(@(a)sqrt(dot(a,a,2)),obj.DLC.tDXYZbp,UniformOutput=false);
			tXYZ_speed=cell2mat(tXYZ_speed(:)');
			tANG_speed=abs(obj.DLC.tDANGbp);
			tDOT_speed=abs(obj.DLC.tDDOTbp);
			obj.DLC_speed.tDXYZbp_magnitude=tXYZ_speed;
			obj.DLC_speed.tDANGbp_magnitude=tANG_speed;
			obj.DLC_speed.tDDOTbp_magnitude=tDOT_speed;
		end
		function [tXYZ_accel,tANG_accel,tDOT_accel]=make_DLC_accel(obj)
			tXYZ_accel=cellfun(@(a)sqrt(dot(a,a,2)),obj.DLC.tD2XYZbp,UniformOutput=false);
			tXYZ_accel=cell2mat(tXYZ_accel(:)');
			tANG_accel=abs(obj.DLC.tD2ANGbp);
			tDOT_accel=abs(obj.DLC.tD2DOTbp);
			obj.DLC_accel.tD2XYZbp_magnitude=tXYZ_accel;
			obj.DLC_accel.tD2ANGbp_magnitude=tANG_accel;
			obj.DLC_accel.tD2DOTbp_magnitude=tDOT_accel;
		end
		function SP_makeBins(obj,options)
			arguments
				obj
				options.phibincount=5
				options.taubincount=6
				options.Xbincount=1
				options.Ybincount=1
				options.Zbincount=5
				options.Xrange_in_mm=[0,0]
				options.Yrange_in_mm=[0,0]
				options.Zrange_in_mm=[-5,5]
			end
			options=namedargs2cell(options);
			[edges_cellVector,bincenter_cellVector]=makeSPBins(options{:});
			obj.SP_edges_cellVector=edges_cellVector;
			obj.SP_bincenter_cellVector=bincenter_cellVector;
		end
		%%
		function SP_Discretize(obj,params2discretize_stringArray4jointDiscretization,options)
			arguments
				obj
			end
			arguments (Repeating)
				params2discretize_stringArray4jointDiscretization (1,:) string {mustBeMember(params2discretize_stringArray4jointDiscretization, ...
					["theta","phi","tau","X","Y","Z"])}
			end
			arguments
				options.Output (1,1) string {mustBeMember(options.Output, ...
					["CellArray","NumericArray_concatenateOutputIn2ndDimension"])}="CellArray"
				options.EncodeBidirectionalPhi logical =true
				options.EncodeBidirectional logical =true
			end
			names={};
			for i=1:length(params2discretize_stringArray4jointDiscretization)
				names{i}=join(params2discretize_stringArray4jointDiscretization{i}(:)',"_",2);
			end

			if options.EncodeBidirectionalPhi
				obj.EncodeBidirectionalPhi=true;
				[newtpk,newedges]=make_tpk_and_edges_alongwith_bidirectional_phi();
				tpk=newtpk;
				edges=newedges;
				obj.tpk_bidirectional=tpk;
			else
				obj.EncodeBidirectionalPhi=false;
				tpk=obj.tpk;
				edges=obj.SP_edges_cellVector;
			end
			if options.EncodeBidirectional
				motionlist=["tau","X","Y","Z"];
				for m=1:length(motionlist)
					[newtpk,newedges]=make_tpk_and_edges_alongwith_bidirectional(motionlist(m),edges);
					tpk=newtpk;
					edges=newedges;
					obj.tpk_bidirectional=tpk;
				end
			end
			[Discretized_t1k_cellVector, DiscretizingFunc_cellVector ]...
				= discretizeSpecifiedSPparams(tpk, edges, params2discretize_stringArray4jointDiscretization{:});
			Func=CellWithStruct2StructWithCell(DiscretizingFunc_cellVector);
			if options.Output=="NumericArray_concatenateOutputIn2ndDimension"
				Discretized_t1k_cellVector=cell2mat(Discretized_t1k_cellVector);
			end
			obj.SP_Discretized_names=names;
			obj.SP_Discretized_cellVector=Discretized_t1k_cellVector;
			obj.SP_Discretized_Discretizer=Func.discretize;
			obj.SP_Discretized_Undiscretizer=Func.undiscretize;
			obj.SP_Discretized_edges=Func.edges;
			obj.SP_Discretized_bincount=Func.bincount;
			obj.SP_Discretized_params2encode=[params2discretize_stringArray4jointDiscretization{:}];
			obj.SP_Discretized_params2encode_cell=params2discretize_stringArray4jointDiscretization;
			obj.SP_Discretized_Undiscretized={};
			for i=1:numel(obj.SP_Discretized_cellVector)
				f=obj.SP_Discretized_Undiscretizer{i};
				obj.SP_Discretized_Undiscretized{i}=f(obj.SP_Discretized_cellVector{i});
			end
			% save original data and the corresponding names
			obj.SP_Discretized_orgData={};
			obj.SP_Discretized_Undiscretized_names={};
			for i=1:numel(params2discretize_stringArray4jointDiscretization)
				IND=makeIND(tpk,2,obj.ind_SPparam(params2discretize_stringArray4jointDiscretization{i}));
				obj.SP_Discretized_orgData{i}=tpk(IND{:});
				obj.SP_Discretized_Undiscretized_names{i}=params2discretize_stringArray4jointDiscretization{i}(:)';
			end

			function [newtpk,newedges]=make_tpk_and_edges_alongwith_bidirectional_phi()
				IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam("phi");
				currdata=obj.tpk(IND{:});
				orgdata=currdata;
				diffdata=diff(currdata,1,1);
				IND=repmat({':'},1,ndims(obj.tpk));IND{1}=size(diffdata,1);
				altIND=IND;altIND{1}=size(diffdata,1)+1;
				diffdata(altIND{:})=diffdata(IND{:});		% make the size same as currdata
				currdata(diffdata==0)=NaN;
				currdata(diffdata>0)=currdata(diffdata>0);
				currdata(diffdata<0)=-currdata(diffdata<0);	% ascending --> positive, descending --> negative
				currdata=fillmissing(currdata,"previous",1);
				currdata=fillmissing(currdata,"next",1);
				currdata=replaceAllNaN_in1stdimension_toOriginal(currdata,orgdata);
				
				phiedges=obj.SP_edges_cellVector{obj.ind_SPparam("phi")};
				nphiedges=-phiedges;
				phiedges=[nphiedges(end:-1:2);phiedges];

				% assign new variables
				newtpk=obj.tpk;IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam("phi");
				newtpk(IND{:})=currdata;
				newedges=obj.SP_edges_cellVector;
				newedges{obj.ind_SPparam("phi")}=phiedges;
			end
			function [newtpk,newedges]=make_tpk_and_edges_alongwith_bidirectional(nameOfMotion,edges)
				IND=repmat({':'},1,ndims(obj.tpk_bidirectional));IND{2}=obj.ind_SPparam(nameOfMotion);
				currdata=obj.tpk_bidirectional(IND{:});
				orgdata=currdata;
				maxval=max(currdata,[],"all");
				minval=min(currdata,[],"all");
				diffdata=diff(currdata,1,1);
				IND=repmat({':'},1,ndims(obj.tpk_bidirectional));IND{1}=size(diffdata,1);
				altIND=IND;altIND{1}=size(diffdata,1)+1;
				diffdata(altIND{:})=diffdata(IND{:});		% make the size same as currdata
				currdata(diffdata==0)=NaN;
				currdata(diffdata>0)=currdata(diffdata>0);
				currdata(diffdata<0)=maxval-minval+currdata(diffdata<0);	% ascending --> same as original, descending --> original + maxval-minval (span of the data)
				currdata=fillmissing(currdata,"previous",1);
				currdata=fillmissing(currdata,"next",1);
				currdata=replaceAllNaN_in1stdimension_toOriginal(currdata,orgdata);
				
				phiedges=edges{obj.ind_SPparam(nameOfMotion)};
				curredges=[minval;phiedges(2:end-1);maxval];
				phiedges=[curredges;curredges(2:end)+maxval-minval];
				phiedges(1)=-Inf;phiedges(end)=Inf;

				% assign new variables
				newtpk=obj.tpk_bidirectional;IND=repmat({':'},1,ndims(obj.tpk_bidirectional));IND{2}=obj.ind_SPparam(nameOfMotion);
				newtpk(IND{:})=currdata;
				newedges=edges;
				newedges{obj.ind_SPparam(nameOfMotion)}=phiedges;
			end
			function outputdata=replaceAllNaN_in1stdimension_toOriginal(nanned_data,orginal_data)
				nancheck=all(isnan(nanned_data),1);
				nanlog=repmat(nancheck,size(nanned_data,1),1);
				nanned_data(nanlog)=orginal_data(nanlog);
				outputdata=nanned_data;
			end
		end
		%%
		function plot_SP_Discretized(obj,param2plot,ntl)
			arguments
				obj
				param2plot
				ntl=[]
			end
			ntl=obj.check_and_assign_nested_tiledlayout(ntl);
			ntl.nexttile();
			f=@(a)BlockCat(reshape(a,size(a,1),size(a,2),[]),3,1);

			Undis_names=obj.SP_Discretized_Undiscretized_names{param2plot};
			tpk=obj.tpk;
			orgdata=tpk(:,obj.ind_SPparam(Undis_names),:,:);
			orgdata=f(orgdata);
			plot(orgdata,DisplayName="Original data");hold on;
			undiscretized=obj.SP_Discretized_Undiscretized{param2plot};
			plot(f(undiscretized),DisplayName="Undiscretized data");
			legend();

			ntl.nexttile();
			data=obj.SP_Discretized_cellVector{param2plot};
			data=f(data);
			plot(data)
			% imagesc(data')
			ntl.change_settings(LinkAxes="x");
			% 
			% currmotions=obj.SP_Discretized_params2encode_cell{param2plot};
			% IND=repmat({':'},1,ndims(obj.tpk));IND{2}=obj.ind_SPparam(currmotions);
			% orgdata=obj.tpk(IND{:});
			% orgdata=f(orgdata);
			% plot(orgdata,DisplayName="Original data");hold on;
			% undiscretized=obj.SP_Discretized_Undiscretized{param2plot};
			% % try
			% % 	func=obj.SP_Discretized_Undiscretizer{param2plot};
			% % 	undiscretized=func(obj.SP_Discretized_cellVector{param2plot});
			% % catch
			% % 	func=obj.SP_Discretized_Undiscretizer{param2plot};
			% % 	undiscretized=func(obj.SP_Discretized_cellVector);
			% % end
			% plot(f(undiscretized),DisplayName="Undiscretized data");
			% legend();
			% 
			% ntl.nexttile();
			% data=obj.SP_Discretized_cellVector{param2plot};
			% data=f(data);
			% plot(data)
			% % imagesc(data')
			% ntl.change_settings(LinkAxes="x");
		end
		%%
		function DLC_makeBins(obj,options)
			% === obj properties used as input ===
			% obj.DLC
			% obj.invDLC
			% obj.rec3d
			% obj.DLC.triError

			% === obj properties written to (modified or set) ===
			% obj.DLC_inv_or_normal
			% obj.DLC_edges_type_cellVector
			% obj.DLC_bincenter_type_cellVector
			% obj.DLC_edges_speed_cellVector
			% obj.DLC_bincenter_speed_cellVector
			% obj.DLC_edges_speed_cellVector_names
			% obj.DLC_edges_accel_cellVector
			% obj.DLC_bincenter_accel_cellVector
			% obj.DLC_edges_accel_cellVector_names
			% obj.DLC_speed
			% obj.DLC_accel

			% === methods called inside ===
			% DLC_NaNfilter_by_triError
			% createBinEdges
			% selectiveDLC_NaNfilter_by_triError

			arguments
				obj
				options.triError_thresh_in_mm	=4	% To use for calculation of bins
				options.postureFormat (1,1) string {mustBeMember(options.postureFormat, ["normalPosture","PostureInversedBySP"])}="normalPosture"
				options.DLC_edge_percent_thresh_For_cutoff=5;
				options.bincount=5;
			end
			triError_thresh_in_mm=options.triError_thresh_in_mm;
			postureFormat=options.postureFormat;
			DLC_edge_percent_thresh_For_cutoff=options.DLC_edge_percent_thresh_For_cutoff;
			bincount=options.bincount;
			obj.DLC_Discretized_bincount=bincount;
			obj.DLC_Discretized_edge_percent_thresh_For_cutoff=options.DLC_edge_percent_thresh_For_cutoff;

			switch postureFormat
				case "normalPosture"
					NaNed_DLC_allcell=DLC_NaNfilter_by_triError(obj.DLC, obj.rec3d, triError_thresh_in_mm,Output="AllNumeric_concatenatedIn2ndDimension");
				case "PostureInversedBySP"
					NaNed_DLC_allcell=DLC_NaNfilter_by_triError(obj.invDLC, obj.rec3d, triError_thresh_in_mm,Output="AllNumeric_concatenatedIn2ndDimension");
			end
			obj.DLC_inv_or_normal=options.postureFormat;
			obj.DLC_edges_type_cellVector=[];
			flds=string(fields(NaNed_DLC_allcell));
			for k=flds(:)'
				currdata=NaNed_DLC_allcell.(k);
				IND=repmat({':'},1,ndims(currdata));
				curredges=cell(1,size(currdata,2));
				currbincenter=cell(1,size(currdata,2));
				for i=1:size(currdata,2)
					% fprintf("%s, %d\n",k,i)
					IND{2}=i;
					currslice=currdata(IND{:});
					currslice=currslice(~isnan(currslice));	% remove NaN values
					[curredges{i},currbincenter{i}]=createBinEdges(currslice,bincount,"Uniform",DLC_edge_percent_thresh_For_cutoff);
				end
				obj.DLC_edges_type_cellVector.(k)=curredges;
				obj.DLC_bincenter_type_cellVector.(k)=currbincenter;
			end

			NaNed_speed=makeEdgesForSpeedMagnitude();
			flds=string(fields(NaNed_speed));
			for k=flds(:)'
				curredges=[];
				currbincenter=[];
				data=NaNed_speed.(k);
				for i=1:size(data,2)
					IND=repmat({':'},1,ndims(data));IND{2}=i;
					currdata=data(IND{:});
					currdata=currdata(~isnan(currdata));
					[curredges{i},currbincenter{i}] = createBinEdges(currdata, bincount, "Uniform", DLC_edge_percent_thresh_For_cutoff);
				end
				obj.DLC_edges_speed_cellVector.(k)=curredges;
				obj.DLC_bincenter_speed_cellVector.(k)=currbincenter;
			end
			obj.DLC_edges_speed_cellVector_names=flds;
			
			NaNed_accel=makeEdgesForAccelMagnitude();
			flds=string(fields(NaNed_accel));
			for k=flds(:)'
				curredges=[];
				currbincenter=[];
				data=NaNed_accel.(k);
				for i=1:size(data,2)
					IND=repmat({':'},1,ndims(data));IND{2}=i;
					currdata=data(IND{:});
					currdata=currdata(~isnan(currdata));
					[curredges{i},currbincenter{i}] = createBinEdges(currdata, bincount, "Uniform", DLC_edge_percent_thresh_For_cutoff);
				end
				obj.DLC_edges_accel_cellVector.(k)=curredges;
				obj.DLC_bincenter_accel_cellVector.(k)=currbincenter;
			end
			obj.DLC_edges_accel_cellVector_names=flds;

			% make eigenposture bins
			if ~isempty(obj.eigen_DLC)
				curredges={};
				for i=1:size(obj.eigen_DLC,2)
					IND=makeIND(obj.eigen_DLC,2,i);
					currdata=obj.eigen_DLC(IND{:});
					[curredges{i},~]=createBinEdges(currdata, bincount, "Uniform", DLC_edge_percent_thresh_For_cutoff);
				end
				obj.eigen_DLC_edges=curredges;
			end

			% generate speed and replace outliers by nan
			function NaNed_speed=makeEdgesForSpeedMagnitude()
				NaNed_speed=[];
				currfields=["tDXYZbp","tDDOTbp","tDANGbp"];
				obj.DLC_speed=[];
				for j=currfields(:)'
					switch postureFormat
						case "normalPosture"
							mdata=obj.DLC;
						case "PostureInversedBySP"
							mdata=obj.invDLC;
					end
					data=mdata.(j);
					if iscell(data)
						speed=cellfun(@(a)sqrt(dot(a,a,2)),data,UniformOutput=false);
						speed=cell2mat(speed);
					else
						speed=abs(data);	% assume they are scalar variables
					end
					% replace by nan
					newdata=selectiveDLC_NaNfilter_by_triError(speed,obj.DLC.triError, obj.rec3d, triError_thresh_in_mm, ...
						Output="AllNumeric_concatenatedIn2ndDimension",XYZOrScalar="Scalar");
					newfield=append(j,"_magnitude");
					NaNed_speed.(newfield)=newdata;
					obj.DLC_speed.(newfield)=speed;
				end
			end
			% generate acceleration and replace outliers by nan
			function NaNed_speed=makeEdgesForAccelMagnitude()
				NaNed_speed=[];
				currfields=["tD2XYZbp","tD2DOTbp","tD2ANGbp"];
				obj.DLC_accel=[];
				for j=currfields(:)'
					switch postureFormat
						case "normalPosture"
							mdata=obj.DLC;
						case "PostureInversedBySP"
							mdata=obj.invDLC;
					end
					data=mdata.(j);
					if iscell(data)
						accel=cellfun(@(a)sqrt(dot(a,a,2)),data,UniformOutput=false);
						accel=cell2mat(accel);
					else
						accel=abs(data);	% assume they are scalar variables
					end
					% replace by nan
					newdata=selectiveDLC_NaNfilter_by_triError(accel,obj.DLC.triError, obj.rec3d, triError_thresh_in_mm, ...
						Output="AllNumeric_concatenatedIn2ndDimension",XYZOrScalar="Scalar");
					newfield=append(j,"_magnitude");
					NaNed_speed.(newfield)=newdata;
					obj.DLC_accel.(newfield)=accel;
				end
			end
		end
	end
	properties
		DLC_Discretized_cellVector
		DLC_Discretized_edges_cellVector
		DLC_Discretized_names
		DLC_Discretized_params2encode
		DLC_Discretized_Undiscretizer
		DLC_Discretized_Undiscretized
		DLC_Discretized_orgData
		DLC_Discretized_Undiscretized_names
		
	end
	methods
		%%
		function DLC_Discretize(obj,params2discretize)
			arguments
				obj
			end
			arguments (Repeating)
				params2discretize (1,1) string {mustBeMember(params2discretize,["tXYZbp","tDXYZbp","tD2XYZbp", ...
					"tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp", ...
					"tDXYZbp_magnitude","tDDOTbp_magnitude","tDANGbp_magnitude", ...
					"tD2XYZbp_magnitude","tD2DOTbp_magnitude","tD2ANGbp_magnitude", ...
					"tDXYZbp_magnitude_mean", ...
					"eigen_DLC"])}
			end
			XYZ=["tXYZbp","tDXYZbp","tD2XYZbp"];
			Ang=["tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"];

			obj.DLC_Discretized_params2encode=[params2discretize{:}];
			obj.DLC_Discretized_cellVector=cell(1,length(params2discretize));
			obj.DLC_Discretized_edges_cellVector={};
			obj.DLC_Discretized_names={};
			% make undiscretizer
			obj.DLC_Discretized_Undiscretizer=cell(1,length(params2discretize));
			obj.DLC_Discretized_Undiscretized=cell(1,length(params2discretize));
			obj.DLC_Discretized_orgData=cell(1,length(params2discretize));
			obj.DLC_Discretized_Undiscretized_names={};

			switch obj.DLC_inv_or_normal
				case "normalPosture"
					DLCdata=obj.DLC_convert(obj.DLC,"AllNumeric");
				case "PostureInversedBySP"
					DLCdata=obj.DLC_convert(obj.invDLC,"AllNumeric");
			end

			for i=1:length(params2discretize)
				currstr=params2discretize{i};
				name_of_current_category=[];
				if isfield(DLCdata,currstr)
					currdata=DLCdata.(currstr);
					curredges=obj.DLC_edges_type_cellVector.(currstr);

					if any(containsSubstring(XYZ,currstr))
						temp=arrayfun(@(a)currstr+"_"+a+["_x","_y","_z"],obj.rec3d.bp_ind,'UniformOutput',false);
						name_of_current_category=cat(2,temp{:});...1 x 42 string array
					elseif any(containsSubstring(Ang,currstr))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif isfield(obj.DLC_speed,currstr)
					currdata=obj.DLC_speed.(currstr);
					curredges=obj.DLC_edges_speed_cellVector.(currstr);

					if any(containsSubstring(currstr,"XYZ"))
						name_of_current_category=currstr+"_"+obj.rec3d.bp_ind;...1 x 14 string array
					elseif any(containsSubstring(currstr,"ANG"))||any(containsSubstring(currstr,"DOT"))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif isfield(obj.DLC_accel,currstr)
					currdata=obj.DLC_accel.(currstr);
					curredges=obj.DLC_edges_accel_cellVector.(currstr);

					if any(containsSubstring(currstr,"XYZ"))
						name_of_current_category=currstr+"_"+obj.rec3d.bp_ind;...1 x 14 string array
					elseif any(containsSubstring(currstr,"ANG"))||any(containsSubstring(currstr,"DOT"))
						name_of_current_category=currstr+"_"+obj.rec3d.jointNames(:)';...1 x 16 string array
					else
						error("No match")
					end
				elseif currstr=="tDXYZbp_magnitude_mean"
					data=obj.DLC_speed.tDXYZbp_magnitude;
					meandata=mean(data,2);
					[curredges,~]=createBinEdges(meandata,obj.DLC_Discretized_bincount,"Uniform",obj.DLC_Discretized_edge_percent_thresh_For_cutoff);
					curredges={curredges};
					currdata=meandata;
					name_of_current_category="tDXYZbp_magnitude_mean";
				elseif currstr=="eigen_DLC"
					if isempty(obj.eigen_DLC)
						error("solve for eigenposture first")
					end
					currdata=obj.eigen_DLC;
					curredges=obj.eigen_DLC_edges;
					name_of_current_category="eigen_posture_"+string(1:size(currdata,2));
				else
					error("Specified param2discretize is not recognized")
				end	
				disc_data=DiscretizeIn2ndDimension(currdata,curredges);
				obj.DLC_Discretized_cellVector{i}=disc_data;
				obj.DLC_Discretized_edges_cellVector{i}=curredges;
				obj.DLC_Discretized_names{i}=cellfun(@string,cellstr(name_of_current_category));
				undiscretized={};
				undiscretizer={};
				for k=1:numel(disc_data)
					undiscretizer{k}=@(discretized)Undiscretize(discretized,curredges{k});
					undiscretized{k}=Undiscretize(disc_data{k},curredges{k});
				end
				obj.DLC_Discretized_Undiscretizer{i}=undiscretizer(:)';
				obj.DLC_Discretized_Undiscretized{i}=undiscretized(:)';
				currsize=num2cell(size(currdata));
				obj.DLC_Discretized_orgData{i}=mat2cell(currdata,size(currdata,1),ones(1,size(currdata,2)),currsize{3:end});
				% make names for undiscretized data
				obj.DLC_Discretized_Undiscretized_names{i}=obj.DLC_Discretized_names{i};
			end
			function disc_data=DiscretizeIn2ndDimension(numeric_data,edges_cell)
				IND=repmat({':'},1,ndims(numeric_data));
				if size(numeric_data,2)~=length(edges_cell);error("Data size and edges length must match");end
				disc_data=[];
				for j=1:size(numeric_data,2)
					IND{2}=j;
					temp=numeric_data(IND{:});
					disc_data{j}=discretize(temp,edges_cell{j});
				end
			end
		end
		%%
		function plot_Discretized_DLC(obj,param2plot,ntl)
			arguments
				obj
				param2plot  (1,1) string
				ntl=[]
			end

			for i=1:length(obj.DLC_Discretized_names)
				ind=obj.DLC_Discretized_names{i}==param2plot;
				if sum(ind)>0
					discretized=obj.DLC_Discretized_cellVector{i}{ind};
					undiscretized=obj.DLC_Discretized_Undiscretized{i}{ind};
					org_data=obj.DLC_Discretized_orgData{i}{ind};
					edges=obj.DLC_Discretized_edges_cellVector{i}{ind};
				end
			end

			% switch obj.DLC_inv_or_normal
			% 	case "normalPosture"
				% 	DLCdata=obj.DLC_convert(obj.DLC,"AllNumeric");
			% 	case "PostureInversedBySP"
				% 	DLCdata=obj.DLC_convert(obj.invDLC,"AllNumeric");
			% end
			% 
			% 
			% currstr=param2plot;
			% if isfield(DLCdata,currstr)
			% 	currdata=DLCdata.(currstr);
			% 	curredges=obj.DLC_edges_type_cellVector.(currstr);
			% elseif isfield(obj.DLC_speed,currstr)
			% 	currdata=obj.DLC_speed.(currstr);
			% 	curredges=obj.DLC_edges_speed_cellVector.(currstr);
			% elseif isfield(obj.DLC_accel,currstr)
			% 	currdata=obj.DLC_accel.(currstr);
			% 	curredges=obj.DLC_edges_accel_cellVector.(currstr);
			% else
			% 	error("No field detected")
			% end
			% discretized=obj.DLC_Discretized_cellVector{containsSubstring(obj.DLC_Discretized_params2encode,currstr)}{ind2plot};
			% edges=obj.DLC_Discretized_edges_cellVector{containsSubstring(obj.DLC_Discretized_params2encode,currstr)}{ind2plot};
			% undiscretized=Undiscretize(discretized,edges);
			% IND=createIND(currdata,2,ind2plot);
			% org_data=currdata(IND{:});
			ntl=obj.check_and_assign_nested_tiledlayout(ntl);
			ntl.nexttile();
			plot(flattenData(org_data));hold on;
			plot(flattenData(undiscretized));
			yline(edges(~isinf(edges)));
			legend(["Original data","Undiscretized data","Edges"])
			ntl.nexttile();
			plot(flattenData(discretized));
			ntl.change_settings(LinkAxes="x");
			function data=flattenData(data)
				data=reshape(permute(data,[2,1,3:ndims(data)]),size(data,2),[])';
			end
		end
		function display_discretized_DLC(obj,param2plot,ind2plot,ntl)
			arguments
				obj
				param2plot  (1,1) string {mustBeMember(param2plot,["tXYZbp","tDXYZbp","tD2XYZbp", ...
					"tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp", ...
					"tDXYZbp_magnitude","tDDOTbp_magnitude","tDANGbp_magnitude", ...
					"tD2XYZbp_magnitude","tD2DOTbp_magnitude","tD2ANGbp_magnitude"])}
				ind2plot
				ntl=[]
			end
			switch obj.DLC_inv_or_normal
				case "normalPosture"
					DLCdata=obj.DLC_convert(obj.DLC,"AllNumeric");
				case "PostureInversedBySP"
					DLCdata=obj.DLC_convert(obj.invDLC,"AllNumeric");
			end
			currstr=param2plot;
			if isfield(DLCdata,currstr)
				currdata=DLCdata.(currstr);
				curredges=obj.DLC_edges_type_cellVector.(currstr);
			elseif isfield(obj.DLC_speed,currstr)
				currdata=obj.DLC_speed.(currstr);
				curredges=obj.DLC_edges_speed_cellVector.(currstr);
			elseif isfield(obj.DLC_accel,currstr)
				currdata=obj.DLC_accel.(currstr);
				curredges=obj.DLC_edges_accel_cellVector.(currstr);
			else
				error("No field detected")
			end
			discretized=obj.DLC_Discretized_cellVector{containsSubstring(obj.DLC_Discretized_params2encode,currstr)}{ind2plot};
			edges=obj.DLC_Discretized_edges_cellVector{containsSubstring(obj.DLC_Discretized_params2encode,currstr)}{ind2plot};
			undiscretized=Undiscretize(discretized,edges);
			IND=createIND(currdata,2,ind2plot);
			org_data=currdata(IND{:});
			ntl=obj.check_and_assign_nested_tiledlayout(ntl);
			ntl.nexttile();
			plot(flattenData(org_data));hold on;
			plot(flattenData(undiscretized));
			yline(edges(~isinf(edges)));
			legend(["Original data","Undiscretized data","Edges"])
			ntl.nexttile();
			plot(flattenData(discretized));
			ntl.change_settings(LinkAxes="x");
			function data=flattenData(data)
				data=reshape(permute(data,[2,1,3:ndims(data)]),size(data,2),[])';
			end
		end
	end
	





%%
	methods (Static)
		function mask=generate_mask4conv_in_plotCellDensity(mask_radius_in_um,mask_shape,padding_in_mm)
			arguments
				mask_radius_in_um
				mask_shape (1,1) string {mustBeMember(mask_shape, ["Gaussian","Ellipse"])}
				padding_in_mm=0
			end
			one_cm=1.77;		% length of one centimeter in powerpoint
			width=18.89;		% length of image width in powerpoint
			height=15.02;		% length of image height in powerpoint
			pixwidth=911;
			pixheight=904;
			x_pix_per_mm=pixwidth*(one_cm/width);
			y_pix_per_mm=pixheight*(one_cm/height);
			mask=generate_mask4conv(x_pix_per_mm,y_pix_per_mm,mask_radius_in_um,mask_shape,padding_in_mm);
		end
	end
%% check data quality
	methods
		function DataQualityCheck_plotRandCellActivity(obj,options)
			arguments
				obj
				options.Seconds2plot=[]
				options.NumberOfCells2plot
				options.skew_thresh
				options.offsetParam=[]
				options.titleFontSize=[]
				options.showTrialBoundary=true
				options.labelFontSize=[]
			end
			tck=obj.tck;
			Skew=obj.Skew;
			IND=repmat({':'},1,ndims(tck));IND{2}=Skew>=options.skew_thresh;
			valtck=tck(IND{:});
			ind=randsample(size(valtck,2),options.NumberOfCells2plot);
			IND{2}=ind(:)';
			randtck=valtck(IND{:});
			if ~isempty(options.Seconds2plot)
				randtck=reshape(randtck,size(randtck,1),size(randtck,2),[]);
				randtck=randtck(:,:,1:round(options.Seconds2plot*10/size(randtck,1)));
			end
			plotCellActivity(randtck, options.offsetParam,options.showTrialBoundary);
			titletext=sprintf("%s, Randomly selected neurons with skewness>=%g",obj.info.NoAltplane.nameFull_p(obj.plane),options.skew_thresh);
			if ~isempty(options.titleFontSize)
				title(titletext,FontSize=options.titleFontSize);
			else
				title(titletext,options.skew_thresh);
			end
			if ~isempty(options.labelFontSize)
				fig=gcf;
				fig.Children.XLabel.FontSize=options.labelFontSize;
				fig.Children.YLabel.FontSize=options.labelFontSize;
			end
		end
		function DataQualityCheck_plotDLCdata(obj,opts,options)
			arguments
				obj
				opts.Data2plot_strarr (1,:) string {mustBeMember(opts.Data2plot_strarr,["tXYZbp","tDXYZbp","tD2XYZbp","tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"])}
				opts.normalize=false
				opts.labelFontSize=[]
				options.offset_prctile=[]
				options.offsetParam=[]
				options.showTrialBoundary=true
				options.labels=[]
				options.labelsize=[]
			end
			DLCdata=obj.DLC_convert(obj.DLC,"AllNumeric");
			dataarray=[];
			namearray=[];
			Ylabel=[];
			Title=[];
			for i=1:length(opts.Data2plot_strarr)
				k=opts.Data2plot_strarr(i);
				Title{i}=k;
				dataarray{i}=DLCdata.(k);
				if opts.normalize
					dataarray{i}=zscore(dataarray{i},0,[1,3:ndims(dataarray{i})]);
				end
				if contains(k,"XYZ")
					namearray{i}=generateXYZnames();
					Ylabel{i}="BodyParts";
				elseif contains(k,"ANG")
					namearray{i}=generateAngnames();
					Ylabel{i}="Joint angles";
				else
					error("No string match")
				end
			end
			obj.plotTile_allPlane(@plotfunc,dataarray,namearray,Ylabel,Title);
			fig=gcf;
			fig.Children.Title.String=obj.info.NoAltplane.nameFull_p(obj.plane);
			if ~isempty(opts.labelFontSize)
				for i=1:numel(fig.Children.Children)
					fig.Children.Children(i).YLabel.FontSize=opts.labelFontSize;
					fig.Children.Children(i).XLabel.FontSize=opts.labelFontSize;
					fig.Children.Children(i).Title.FontSize=opts.labelFontSize;
					fig.Children.Children(i).Title.FontWeight="bold";
					fig.Children.Title.FontSize=opts.labelFontSize;
					fig.Children.Title.FontWeight="bold";
				end
			end
			function plotfunc(data,labelarray,Ylabel,Title)
				optns=options;
				optns.labels=labelarray;
				optns=namedargs2cell(optns);
				plotSimilarActivity(data,optns{:});
				ylabel(Ylabel)
				title(Title)
			end
			function names=generateXYZnames()
				names=obj.rec3d.name_ind;
				names=repelem(names,1,3);
				xyz=repmat([" x"," y"," z"],1,length(obj.rec3d.name_ind));
				names=append(names,xyz);
			end
			function names=generateAngnames()
				names=obj.rec3d.Dotname_ind(:,1);
			end
		end
		function DataQualityCheck_plotcompare_triError_and_DLCvalues(obj,opts)
			arguments
				obj
				opts.Data2plot_strarr (1,:) string {mustBeMember(opts.Data2plot_strarr,["tXYZbp","tDXYZbp","tD2XYZbp","tDOTbp","tDDOTbp","tD2DOTbp","tANGbp","tDANGbp","tD2ANGbp"])}
				opts.normalize=false
				opts.labelFontSize=[]
			end
			DLCdata=obj.DLC_convert(obj.DLC,"AllNumeric");
			dataarray=[];
			triArray=[];
			namearray=[];
			Xlabel=[];
			Ylabel=[];
			Title=[];
			for i=1:length(opts.Data2plot_strarr)
				k=opts.Data2plot_strarr(i);
				Title{i}=k;
				data=DLCdata.(k);
				data=zscore(data,0,[1,3:ndims(data)]);
				if contains(k,"XYZ")
					data=generateDistance(data);
					dataarray{i}=data;
					triArray{i}=obj.DLC.triError;
					if ~isequal(size(data), size(obj.DLC.triError))
						error('Size must match');
					end
					dataarray{i}=data(:);
					triArray{i}=obj.DLC.triError(:);
					Xlabel{i}="Distance fom mean coordinates (normalized)";
					Ylabel{i}="Triangulation Error (mm)";
					% namearray{i}=obj.rec3d.name_ind;
					% Ylabel{i}="BodyParts";
				elseif contains(k,"ANG")
					data=abs(data);
					triArray{i}=generateTriError();
					if ~isequal(size(data), size(triArray{i}))
						error('Size must match');
					end
					dataarray{i}=data(:);
					triArray{i}=triArray{i}(:);
					Xlabel{i}="Distance fom mean degrees (normalized)";
					Ylabel{i}="Triangulation Error (mm)";
					% namearray{i}=generateAngnames();
					% Ylabel{i}="Joint angles";
				else
					error("No string match")
				end
			end

			obj.plotTile_allPlane(@plotfunc,dataarray,triArray,Title,Xlabel,Ylabel);
			fig=gcf;
			fig.Children.Title.String=obj.info.NoAltplane.nameFull_p(obj.plane);
			if ~isempty(opts.labelFontSize)
				for i=1:numel(fig.Children.Children)
					fig.Children.Children(i).YLabel.FontSize=opts.labelFontSize;
					fig.Children.Children(i).XLabel.FontSize=opts.labelFontSize;
					fig.Children.Children(i).Title.FontSize=opts.labelFontSize;
					fig.Children.Children(i).Title.FontWeight="bold";
					fig.Children.Title.FontSize=opts.labelFontSize;
					fig.Children.Title.FontWeight="bold";
				end
			end

			function plotfunc(data,trierror,Title,Xlabel,Ylabel)
				scatter(data,trierror,10,[1,1,1]*0.3,"filled");
				xlabel(Xlabel);
				ylabel(Ylabel)
				title(Title);
			end

			function celldata=convertMat2cell(data)
				celldata=[];
				for j=1:size(data,2)
					IND=repmat({':'},1,ndims(data));IND{2}=j;
					celldata{j}=data(IND{:});
				end
			end
			function distar=generateDistance(data)
				distar=[];
				for j=1:size(data,2)/3
					IND=repmat({':'},1,ndims(data));IND{2}=[3*(j-1)+1:3*j];
					temp=data(IND{:});
					dist=sqrt(sum(temp.^2,2));
					distar=cat(2,distar,dist);
				end
			end
			function combinedError=generateTriError()
				combinedError=[];
				for j=1:size(obj.rec3d.Dotname_ind,1)
					currstr=obj.rec3d.Dotname_ind(j,:);
					currind=obj.rec3d.ind_name(currstr);
					IND=repmat({':'},1,ndims(obj.DLC.triError));IND{2}=currind;
					currError=obj.DLC.triError(IND{:});
					currError=max(currError,[],2);
					combinedError=cat(2,combinedError,currError);
				end
			end
		end
%% plot to make sure SP and DLC is in sync for p8
		function plotSPandDLC(obj,options)
			arguments
				obj
				options.ShowLegend=true
			end
			direction=8;
			BodyParts2plot=["RightToe","LeftToe"];
			DLCdata=obj.DLC.tXYZbp;
			DLCdata=DLCdata(obj.rec3d.ind_name(BodyParts2plot));
			Z=[];
			XY=[];
			for i=1:length(DLCdata)
				IND=createIND(DLCdata{i},2,3);
				IND{4}=8;
				Z{i}=squeeze(DLCdata{i}(IND{:}));	% time x trials
			end
			IND=createIND(obj.tpk,2,2);
			IND{4}=8;
			SPdata=60*squeeze(obj.tpk(IND{:}));	% time x trials
			h={};			
			col=["red","blue"];
			for i=1:length(Z)
				h{end+1}=plot(Z{i},"Color",col(i),LineWidth=0.5,DisplayName=BodyParts2plot(i));
				hold on;
			end
			h{end+1}=plot(SPdata,"Color",[1,1,1]*0.2,LineWidth=2,DisplayName="SP phi");
			hold off;
			H=gobjects();
			for i=1:length(h)
				H(i)=h{i}(1);
			end
			if options.ShowLegend
				legend(H);
			end
		end
		%%
		function plotSPandDLC_allPlanes(obj,options)
			arguments
				obj
				options.Funcs_loadData_InputIs_csp_and_plane=[]
			end
			for i=1:length(obj.info.NoAltplane.nameFull_p)
				nexttile;
				csp=copy(obj);
				if ~isempty(options.Funcs_loadData_InputIs_csp_and_plane)
					executeFuncsSequentially(csp,i,options.Funcs_loadData_InputIs_csp_and_plane);
				else
					csp.load_specified(i,"p8_tpkd","p8_DLC");
				end
				tpk=PreprocessSPdata(csp.p8_tpkd);
				csp.tpk=tpk;
				csp.DLC=csp.p8_DLC;
				if i==length(obj.info.NoAltplane.nameFull_p)
					csp.plotSPandDLC(ShowLegend=true);
				else
					csp.plotSPandDLC(ShowLegend=false);
				end
				title(csp.info.NoAltplane.nameFull_p(i));
				i
				drawnow;
			end
			ChangeFigureSettings(LegendLoc="east",LegendSize=15);
			function executeFuncsSequentially(csp,plane,funcs)
				for j=1:length(funcs)
					f=funcs{j};
					[nInputs, inputNames] = getFuncHandleInputs(f);
					checkInputFormat(inputNames);
					args={csp,plane};
					f(args{1:nInputs});
				end
				function checkInputFormat(inputNames)
					orgstr=["csp","plane"];
					for k=1:length(inputNames)
						if inputNames(k)~=orgstr(k)
							error("function handle must have inputs of (csp,plane) or (csp)")
						end
					end
				end
			end
		end
		
	end

	%% cPCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		DimReduct
	end
	methods
		function initialize_DimReduct(obj,data2use,NameOfData)
			obj.DimReduct=[];
			obj.DimReduct.(NameOfData)=data2use;
		end
		function data=get_DimReduct(obj,PropStringArray)
			Props=obj.strArray2cellOfStrings(PropStringArray);
			data=getfield(obj.DimReduct,Props{:});
		end
		function set_DimReduct(obj,PropStringArray,data)
			Props=obj.strArray2cellOfStrings(PropStringArray);
			obj.DimReduct=setfield(obj.DimReduct,Props{:},data);
		end
		function cellOfStrings=strArray2cellOfStrings(obj,strArray)
			temp=cellstr(strArray);
			cellOfStrings=cellfun(@string,temp,UniformOutput=false);
		end
		function strArray=append2lastProp(obj,strArray,str2append)
			strArray(end)=strArray(end)+str2append;
		end
		function [new_tck,pca_var]=applyPCA(obj,options)
			arguments
				obj
				options.Prop2apply (1,:) string
				options.numDim2keep=[];
			end
			tck=obj.get_DimReduct(options.Prop2apply);
			tc=multiBlockCat(tck,[3:ndims(tck)],1);
			if isempty(options.numDim2keep)
				options.numDim2keep=size(tc,2);
			end
			[Axes,~,pca_var]=pca(tc,NumComponents=options.numDim2keep);
			new_tck=pagemtimes(tck,Axes);
			obj.set_DimReduct(obj.append2lastProp(options.Prop2apply,"_pca"),new_tck);
			obj.set_DimReduct(obj.append2lastProp(options.Prop2apply,"_pcaAxes"),Axes);
		end
		function applyNFDA(obj,tck)
			arguments
				obj
				tck
			end
			[coeff,Score,A,B,ScoreA,ScoreB,val]=nFDA(td_c_k,'nFDA',weight);
		end
		function applyPCAandNFDA_to_tckd(obj,options)
			arguments
				obj
				options.numDim2keep=[];
			end
			obj.applyPCA(obj.tck);
		end
	end






	%% Miscellaneous
	methods (Static)
		function Fig=figureTile(rowcol,options)
			arguments
				rowcol.rownum
				rowcol.colnum
				options.TileSpacing	(1,1) string {mustBeMember(options.TileSpacing, ["loose","compact","tight","none"])}="tight"
				options.Padding		(1,1) string {mustBeMember(options.Padding, ["loose","compact","tight"])}="tight"
			end
			nvCell = namedargs2cell(options);
			Fig=figure;
			Fig.Units="normalized";
			Fig.Position=[0,0,1,1];
			try
				t=tiledlayout(rowcol.rownum,rowcol.colnum,nvCell{:});
			catch
				t=tiledlayout('flow',nvCell{:});
			end
		end
		%%
		function plotTile_allPlane(plot_func_handle,Data_cellArray,options)
			arguments
				plot_func_handle
			end
			arguments (Repeating)
				Data_cellArray
			end
			arguments
				options.linkaxes (1,1) string {mustBeMember(options.linkaxes,["x","y","xy",""])}=""
				options.FontSize=[]
				options.Legend (1,:) string =""
				options.LegendLoc (1,1) string {mustBeMember(options.LegendLoc,["north","south","east","west","northeast","northwest","southeast","southwest", ...
					"northoutside","southoutside","eastoutside","westoutside","northeastoutside","northwestoutside","southeastoutside","southwestoutside","best","bestoutside","layout",""])}=""
				options.TileHandle=[]
				options.Exec_nexttile=true
			end
			Data=Data_cellArray;
			if sum(abs(diff(cellfun(@(a)length(a),Data))))~=0
				error('cell arrays must have the same size')
			end
			ax=gobjects(0);
			for i=1:length(Data{1})
				if options.Exec_nexttile
					if isempty(options.TileHandle)
						ax(i)=nexttile;
					else
						ax(i)=nexttile(options.TileHandle);
					end
				end
				currdata=cell(1,length(Data));
				for j=1:length(Data)
					try
						currdata{j}=Data{j}{i};
					catch
						currdata{j}=Data{j}(i);
					end
				end
				plot_func_handle(currdata{:});
				% title(obj.info.NoAltplane.nameFull_p);
			end
			fig=gcf;
			if options.linkaxes~=""
				linkaxes(ax,options.linkaxes);
			end
			if ~isempty(options.FontSize)
				fig.Children.Title.FontSize=options.FontSize;
			end
			if options.Legend~=""
				lgd=legend(options.Legend);
				if options.LegendLoc~=""
					lgd.Layout.Tile=options.LegendLoc;
				end
			end
		end

		%%
		function plotSubPlot_allPlane(plot_func_handle,Data_cellArray,options)
			arguments
				plot_func_handle
			end
			arguments (Repeating)
				Data_cellArray
			end
			arguments
				options.linkaxes (1,1) string {mustBeMember(options.linkaxes,["x","y","xy",""])}=""
				options.rownum_colnum (1,:) =[]
			end
			Data=Data_cellArray;
			if sum(abs(diff(cellfun(@(a)length(a),Data))))~=0
				error('cell arrays must have the same size')
			end
			% obtain rownum and colnum
			if ~isempty(options.rownum_colnum)
				rowcol={options.rownum_colnum(1),options.rownum_colnum(2)};
			else
				b=ceil(sqrt(9/16*length(Data{1})));
				a=ceil(length(Data{1})/b);
				rowcol={b,a};
			end

			ax=gobjects(0);
			for i=1:length(Data{1})
				ax(i)=subplot(rowcol{:},i);
				currdata=cell(1,length(Data));
				for j=1:length(Data)
					try
						currdata{j}=Data{j}{i};
					catch
						currdata{j}=Data{j}(i);
					end
				end
				plot_func_handle(currdata{:});
				% title(obj.info.NoAltplane.nameFull_p);
			end
			if options.linkaxes==""
				linkaxes(ax,options.linkaxes);
			end
		end
		%%
		function Fig=figureTileplot(plot_func_handle,Data_cellArray,rowcol,options)
			arguments
				plot_func_handle
			end
			arguments (Repeating)
				Data_cellArray	(1,:){mustBeMember(Data_cellArray,"cell")}
			end
			arguments
				rowcol.rownum
				rowcol.colnum
				options.TileSpacing	(1,1) string {mustBeMember(options.TileSpacing, ["loose","compact","tight","none"])}="tight"
				options.Padding		(1,1) string {mustBeMember(options.Padding, ["loose","compact","tight"])}="tight"
				options.linkaxes	(1,1) string {mustBeMember(options.linkaxes, ["x","y","xy"])}
			end
			nvCell = namedargs2cell(options);
			Fig=figure;
			Fig.Units="normalized";
			Fig.Position=[0,0,1,1];
			try
				t=tiledlayout(rowcol.rownum,rowcol.colnum,nvCell{:});
			catch
				t=tiledlayout('flow',nvCell{:});
			end
			Data=Data_cellArray;
			if sum(abs(diff(cellfun(@(a)length(a),Data))))~=0
				error('cell arrays must have the same size')
			end
			for i=1:length(Data{1})
				nexttile;
				currdata=cell(1,length(Data));
				for j=1:length(Data)
					currdata{j}=Data{j}{i};
				end
				plot_func_handle(currdata{:});
			end
		end
		
			
	end
end
