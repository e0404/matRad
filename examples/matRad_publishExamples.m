% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad publish examples script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This file will publish all examples


%% clear environment
clc,clear,close all;

% make sure you are in the matRad root directory and have added the
% exampels folder to the matlab path

format     = {'pdf'}; %{'html','pdf'};
matRadRoot = pwd;
folderInfo = dir([matRadRoot filesep 'examples']);

   
for ixFile = 1:size(folderInfo,1)
      
    for ixFormat = 1:numel(format)
        
        %make sure not to use auto saved files (*.asv)
        [pathstr,name,ext] = fileparts(folderInfo(ixFile).name);
        
       if ~folderInfo(ixFile).isdir && ~strcmp(name,mfilename) && strcmp(ext,'.m') ...
           
           save('workspace')
           publish(folderInfo(ixFile).name,format{ixFormat});
           clc,clear,close all
           load('workspace')
           
       end
    end
    
end

delete workspace.mat













