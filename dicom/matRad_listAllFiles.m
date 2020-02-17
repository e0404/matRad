function fileList = matRad_listAllFiles(dirPath,uiInput)
% matRad function to get all fiiles in arbitrary deep subfolders
% 
% call
%   fileList = matRad_listAllFiles(dirPath,uiInput)
%
% input
%   dirPath:        (optional) initial folder to start searching
%   uiInput:        (optional) if userInteraction is wanted.
%                   Use EITHER dirPath or uiInput
%
% output
%   fileList:       Filelist
%
% References
%   -
% Note
%                   MATLAB has an internal recursion limit. If you want to
%                   list a huge amount of files you may have to: 
%                   >> % N = Number of allowed recursions
%                   >> set(0,'RecursionLimit',N)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check whether input argument dirPath is parsed
    if ~exist('dirPath','var')
        if exist('uiInput','var')
            if uiInput == 1
                dirPath = uigetdir;
            end
        else
            dirPath = pwd;
        end
    end
    
%% get all files   
    
    % info of main directory
    mainDirInfo = dir(dirPath);
    % index of subdirectories
    dirIx = [mainDirInfo.isdir];
    % list all files which are not subDirs
    fileList = {mainDirInfo(~dirIx).name}';
    if ~isempty(fileList)
        % using fullfile to take care on file seperators
        fileList = cellfun(@(x) fullfile(dirPath,x),...
                           fileList,'UniformOutput',false);
    end
    % subdirectories
    subDirList = {mainDirInfo(dirIx).name};
    % list items without '.' or '..'
    validIx = ~ismember(subDirList,{'.','..'});
    % loop subDir and recall function itself
    for ixDir = find(validIx)
        % get subDir Path
        nextSubDir = fullfile(dirPath,subDirList{ixDir});
        fileList = [fileList; matRad_listAllFiles(nextSubDir)];
    end

end
