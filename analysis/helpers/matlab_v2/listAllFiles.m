function fileList = listAllFiles(folderPath)
% LISTFILES Return absolute paths of all files under a folder (including subfolders)
%   fileList = LISTFILES(folderPath) returns a string array of full file paths.

    % Validate input
    if ~isfolder(folderPath)
        error('"%s" is not a valid folder.', folderPath);
    end

    % Recursively list everything, then keep only files
    info = dir(fullfile(folderPath, '**', '*'));
    info([info.isdir]) = [];  

    % Build absolute paths and convert to string array
    fileList = string(fullfile({info.folder}, {info.name}))';
end