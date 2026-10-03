% Export the registered ROI masks of the eight cortical regions used by the
% Fig. 6b analysis (BrainRegions_with_RSC.mat) to region_masks_with_RSC.mat.
% Save the file in the results folder that fig6b_regional_lri.py reads
% (analysis/results).

path2functions        = '/path/to/matlab_helpers';
path2allMiceDataDirlist = "/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat";
path2BrainRegions     = "/path/to/data/BrainRegions_with_RSC.mat";
outfile = "/path/to/analysis/results/region_masks_with_RSC.mat";

addpath(genpath(path2functions));

nPlane   = numel(CSP().PlaneName);
maskCell = cell(1, nPlane);
planeNames = strings(1, nPlane);
regNames = [];

for plane = 1:nPlane
    csp = CSP(dir_AllMiceDataFile_used=path2allMiceDataDirlist, ...
              dir_BrainRegions_used=path2BrainRegions);
    csp.load_specified(plane, "regROI");
    logs = csp.Reg_inRegLog_reg;              % cell array, one logical per region
    M = false(numel(logs{1}), numel(logs));
    for r = 1:numel(logs)
        M(:, r) = logical(logs{r}(:));
    end
    maskCell{plane} = M;
    planeNames(plane) = csp.PlaneName(plane);
    if isempty(regNames)
        regNames = convertCharsToStrings(csp.Reg_RegionNames);
    end
    fprintf('plane %2d/%d %-16s nROI=%5d  perRegion=%s\n', ...
        plane, nPlane, planeNames(plane), size(M,1), mat2str(sum(M,1)));
end

regNamesChar = cellstr(regNames);
planeNamesChar = cellstr(planeNames);
save(outfile, "maskCell", "regNamesChar", "planeNamesChar", "-v7");
fprintf('\nwrote %s\n', outfile);
