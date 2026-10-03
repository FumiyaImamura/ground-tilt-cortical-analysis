function fit_theta_phi_ablation(planeIndices, fitFullModel)
% Fit the full encoding model and the tilt-predictor ablation.
%
% The full model contains theta, phi, tau, Z, posture and
% movement.  This reduced model removes only theta and phi while retaining
% tau, Z, posture and movement. Both fits use the same preprocessing.
%
% Examples:
%   fit_theta_phi_ablation(1, true)    % full and reduced models, one session
%   fit_theta_phi_ablation(1:20,false) % all reduced models, resume-safe

    arguments
        planeIndices (1,:) double {mustBeInteger,mustBePositive} = 1:20
        fitFullModel (1,1) logical = true
    end

    path2functions = '/path/to/matlab_helpers';
    path2allMiceDataDirlist = "/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2.mat";
    path2BrainRegions = "/path/to/data/BrainRegions.mat";
    outputDir = "/path/to/encoding/results/theta_phi_ablation_scores";
    if ~isfolder(outputDir); mkdir(outputDir); end
    addpath(genpath(path2functions));

    sigma = 1;
    ridgeparams = 10.^(2:0.5:8);
    groupnum = 10;
    make_Data_single = false;

    template = CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
        dir_BrainRegions_used=path2BrainRegions);
    nPlane = numel(template.PlaneName);
    assert(all(planeIndices <= nPlane), 'Requested plane exceeds available planes.');

    for p = planeIndices
        planeName = string(template.PlaneName4Directory(p));
        outputFile = fullfile(outputDir, planeName + ".mat");
        if isfile(outputFile)
            fileInfo = whos('-file',outputFile);
            vars = string({fileInfo.name});
            if any(vars == "ThetaPhiAblated_score_ridge_cells")
                fprintf('Skipping completed output %d/%d\n',p,nPlane);
                continue;
            end
        end

        fprintf('\n=== theta/phi ablation %d/%d ===\n',p,nPlane);
        planeTimer = tic;
        csp = CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
            dir_BrainRegions_used=path2BrainRegions);
        csp.loadData(p,PeriodCutoff4Derivatives_sec=0.1,triErrorThresh_mm=4);
        csp.loadData_continuous_time(UseCachedData=true);
        csp.setCurrentData_to_continuous();

        % Revised reduced model: retain tau and Z; omit theta and phi.
        csp.SP_makeBins(phibincount=5,taubincount=5,Zbincount=5);
        csp.SP_Discretize("tau","Z", ...
            EncodeBidirectionalPhi=true,EncodeBidirectional=true);
        csp.SP_GaussianHotEncode_Discretized(sigma);
        csp.DLC_makeBins(bincount=5,DLC_edge_percent_thresh_For_cutoff=5, ...
            postureFormat="normalPosture",triError_thresh_in_mm=4);
        csp.DLC_Discretize("tXYZbp","tDXYZbp","tD2XYZbp", ...
            "tANGbp","tDANGbp","tD2ANGbp", ...
            "tDXYZbp_magnitude","tDANGbp_magnitude", ...
            "tD2XYZbp_magnitude","tD2ANGbp_magnitude");
        csp.DLC_GaussianHotEncode_Discretized(sigma);
        csp.combinePredictors("Both",true,removeSmilarParameters=false);

        reduced_predictor_count = size(csp.all_GaussianHot,2);
        reduced_predictor_names = csp.all_GaussianHot_names;
        sample_count = size(csp.all_GaussianHot,1);
        roi_count = size(csp.cellac_continuous,2);
        fprintf('Reduced design: %d samples x %d predictors; %d ROIs\n', ...
            sample_count,reduced_predictor_count,roi_count);

        csp.split_Data_to_xvalidation_groups_continuous( ...
            Predictor=csp.all_GaussianHot, ...
            Responser=csp.cellac_continuous, ...
            groupnum=groupnum,TimeChunkSec=5);
        csp.RidgeDecode(Predictor2Use="Predictor",Responser2Use="Responser", ...
            make_Data_single=make_Data_single,ridgeparams=ridgeparams,UseGPU=false);
        ThetaPhiAblated_score_ridge_cells = struct('R', ...
            csp.xval_assess_R2_ridge(ReferenceResponser="Responser"));

        % RidgeDecode retains the fold-level predictor/response/reconstruction
        % arrays in the CSP object.  Release them before building the matched
        % full model so large sessions do not hold both fits in memory.
        clear csp

        RefitFull_score_ridge_cells = [];
        full_predictor_count = NaN;
        full_predictor_names = strings(0,1);

        if fitFullModel
            fprintf('Fitting the full model...\n');
            cspFull = CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
                dir_BrainRegions_used=path2BrainRegions);
            cspFull.loadData(p,PeriodCutoff4Derivatives_sec=0.1,triErrorThresh_mm=4);
            cspFull.loadData_continuous_time(UseCachedData=true);
            cspFull.setCurrentData_to_continuous();
            cspFull.SP_makeBins(phibincount=5,taubincount=5,Zbincount=5);
            cspFull.SP_Discretize(["theta","phi"],"tau","Z", ...
                EncodeBidirectionalPhi=true,EncodeBidirectional=true);
            cspFull.SP_GaussianHotEncode_Discretized(sigma);
            cspFull.DLC_makeBins(bincount=5,DLC_edge_percent_thresh_For_cutoff=5, ...
                postureFormat="normalPosture",triError_thresh_in_mm=4);
            cspFull.DLC_Discretize("tXYZbp","tDXYZbp","tD2XYZbp", ...
                "tANGbp","tDANGbp","tD2ANGbp", ...
                "tDXYZbp_magnitude","tDANGbp_magnitude", ...
                "tD2XYZbp_magnitude","tD2ANGbp_magnitude");
            cspFull.DLC_GaussianHotEncode_Discretized(sigma);
            cspFull.combinePredictors("Both",true,removeSmilarParameters=false);
            full_predictor_count = size(cspFull.all_GaussianHot,2);
            full_predictor_names = cspFull.all_GaussianHot_names;
            assert(size(cspFull.all_GaussianHot,1)==sample_count, ...
                'Full and reduced models have different sample counts.');
            assert(size(cspFull.cellac_continuous,2)==roi_count, ...
                'Full and reduced models have different ROI counts.');

            cspFull.split_Data_to_xvalidation_groups_continuous( ...
                Predictor=cspFull.all_GaussianHot, ...
                Responser=cspFull.cellac_continuous, ...
                groupnum=groupnum,TimeChunkSec=5);
            cspFull.RidgeDecode(Predictor2Use="Predictor",Responser2Use="Responser", ...
                make_Data_single=make_Data_single,ridgeparams=ridgeparams,UseGPU=false);
            RefitFull_score_ridge_cells = struct('R', ...
                cspFull.xval_assess_R2_ridge(ReferenceResponser="Responser"));

            clear cspFull
        end

        model_specification = struct( ...
            'full_environmental_variables',{{'theta','phi','tau','Z'}}, ...
            'reduced_environmental_variables',{{'tau','Z'}}, ...
            'removed_environmental_variables',{{'theta','phi'}}, ...
            'gaussian_sigma_bins',sigma,'cross_validation_folds',groupnum, ...
            'time_chunk_seconds',5,'ridge_parameters',ridgeparams);

        tmpFile = outputFile + ".partial";
        save(tmpFile,'ThetaPhiAblated_score_ridge_cells', ...
            'RefitFull_score_ridge_cells','ridgeparams','model_specification', ...
            'reduced_predictor_count','reduced_predictor_names', ...
            'full_predictor_count','full_predictor_names','sample_count','roi_count','-v7.3');
        movefile(tmpFile,outputFile,'f');
        fprintf('Saved %s [%.1f min]\n',outputFile,toc(planeTimer)/60);
        clear ThetaPhiAblated_score_ridge_cells RefitFull_score_ridge_cells
    end
end
