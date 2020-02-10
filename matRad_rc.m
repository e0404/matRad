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


% set log level accordingly if you do _not_ want to do unit testing
if ~exist('unitTestBool','var') || ~unitTestBool

    param.calcDoseDirect = false;
    param.subIx          = [];
    param.logLevel       = 1;
    
% set log level accordingly if want to do unit testing
else
    
    param.calcDoseDirect = false;
    param.subIx          = [];
    param.logLevel       = 3;
    
end