function [tmpPath,status] = helper_temporaryFolder(folderName,clearIfExists)
%helper_temporaryFolder Creates a temporary folder for test data in the
%  users temporary systemdirectory. 

matRad_cfg = MatRad_Config.instance();
if matRad_cfg.isOctave
    confirm_recursive_rmdir(false,"local");
end

if nargin < 2
    clearIfExists = true;
end

tmpPath = fullfile(tempdir(),folderName);

if isfolder(tmpPath) && clearIfExists    
    status = rmdir(tmpPath,'s');    
    if status
        status = mkdir(tmpPath);
    else
        status = double(isfolder(tmpPath));
    end
elseif ~isfolder(tmpPath)
    status = mkdir(tmpPath);
else
    status = 1;
end

if ~status
    tmpPath = pwd();
    warning('temporary directory invalid - using wokring directory!');
end


end

