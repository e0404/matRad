function matRad_cfg = matRad_rc(clearWindow)
% matRad rc script
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    clearWindow = false;
end

thisFolder = fileparts(mfilename("fullpath"));

sameMatRad = false;
% Initialize matRad
try %Check if matRad is already on the path and if so, if it is the same installation
    tmp_cfg = MatRad_Config.instance();        
    if ~strcmp(tmp_cfg.matRadRoot,thisFolder)
        tmp_cfg.dispWarning('Called matRad_rc in folder %s but matRad initialized in folder %s!\n Removing old matRad from path and using new installation!',fileparts(mfilename("fullpath")),tmp_cfg.matRadRoot);
        if ~isdeployed
            rmpath(genpath(tmp_cfg.matRadRoot));
        end
        clear tmp_cfg MatRad_Config;
        tmp_cfg = MatRad_Config.instance();
        sameMatRad = false;
    else
        sameMatRad = true;
    end
catch %If matRad is not on the path, initialize it freshly and add the sourcefolder containing MatRad_Config to the path
    matRadRoot = thisFolder;
    matRadSrcRoot = fullfile(thisFolder,'matRad');
    if ~isdeployed
        addpath(matRadRoot,matRadSrcRoot);
    end
    tmp_cfg = MatRad_Config.instance();
    sameMatRad = false;
end

% clear workspace and command prompt, close all figures
if clearWindow
    clc;
    close all;
end

% Version Display
vString = matRad_version();

if ~sameMatRad || clearWindow
    %Version Info
    versionStr = sprintf('You are running matRad %s with %s %s',vString,tmp_cfg.env,tmp_cfg.envVersion);
    infoStr = matRad_info();

    message = [versionStr newline infoStr newline];

    tmp_cfg.dispInfo(message);
end

%Return MatRad_Config instance
if nargout > 0
    matRad_cfg = tmp_cfg;
else
    assignin('base','matRad_cfg',tmp_cfg);
end