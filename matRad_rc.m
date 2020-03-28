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

% set search path
addpath(genpath(pwd));

%clear command window and close all figures
clc;
close all;

% clear workspace and command prompt, close all figures
[env,envver] = matRad_getEnvironment();
vString = matRad_version();

fprintf('You are running matRad %s with %s %s\n',vString,env,envver);

switch env
    case 'MATLAB'
        clearvars -except unitTestBool
    case 'OCTAVE'
        clear -x unitTestBool
end

global matRad_cfg;
matRad_cfg = MatRad_Config.instance();

% set log level according to unit testing
if exist('unitTestBool','var') && unitTestBool
    matRad_cfg.logLevel   = 1;
    matRad_cfg.propStf.defaultLongitudinalSpotSpacing = 20;
    matRad_cfg.propDoseCalc.defaultResolution = struct('x',5,'y',6,'z',7); %[mm]
    matRad_cfg.propDoseCalc.defaultGeometricCutOff = 20;
    matRad_cfg.propDoseCalc.defaultLateralCutOff = 0.8;
    matRad_cfg.propOpt.defaultMaxIter = 10;
    matRad_cfg.propMC.ompMC_defaultHistories = 100;
    matRad_cfg.propMC.MCsquare_defaultHistories = 100;
    matRad_cfg.propMC.direct_defaultHistories = 100;
end