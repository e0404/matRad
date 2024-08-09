function [subFolder] = matRad_getSubfolders(folderName)
% matRad_getSubfolders: 
% Recursively loops into the folder and subfolders and outputs the list of
% folders and subfolder names
% call
%   [subFolder] = matRad_getSubfolders(folderName)
%    input:
%       folderName
%   output:
%       subFolders
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the
% help edit
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subFolder = {folderName};
currSubFolders = dir(folderName);

% Filter normal files
currSubFolders = currSubFolders([currSubFolders.isdir]);
currSubFolders = currSubFolders(~ismember({currSubFolders(:).name}, {'.', '..'}));
    
for subIdx=1:numel(currSubFolders)
    subFolder = [subFolder, matRad_getSubfolders(fullfile(currSubFolders(subIdx).folder, currSubFolders(subIdx).name))];
end

end