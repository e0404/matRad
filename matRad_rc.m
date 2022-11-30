% matRad rc script
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

% Initialize matRad
matRad_cfg = MatRad_Config.instance();
if ~strcmp(matRad_cfg.matRadRoot,fileparts(mfilename("fullpath")))
    matRad_cfg.dispWarning('Called matRad_rc in folder %s but matRad initialized in folder %s!\n Removing old matRad from path and using new installation!',fileparts(mfilename("fullpath")),matRad_cfg.matRadRoot); 
    rmpath(genpath(matRad_cfg.matRadRoot));
    clear matRad_cfg MatRad_Config;
    matRad_cfg = MatRad_Config.instance();
end

%clear command window and close all figures
clc;
close all;

% clear workspace and command prompt, close all figures
[env,envver] = matRad_getEnvironment();
vString = matRad_version();

matRad_cfg.dispInfo('You are running matRad %s with %s %s\n',vString,env,envver);
clear env envver vString;


