function matRad_runTests(folder)
%% matRad_runTests.m
% This function runs the test suite for the matRad package.
%
% Usage:
%   - Run the script to execute all the tests in the test suite.
%
% Notes:
%   - Make sure to have MOxUnit installed and added to the path.
%
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

if nargin == 0
    folder = 'test';
end

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispInfo('Setting default properties for testing and starting test suite!\n');
matRad_cfg.setDefaultPropertiesForTesting();
matRad_cfg.logLevel = 1;

back = cd(matRad_cfg.matRadRoot);
moxunit_runtests(folder,'-recursive');

matRad_cfg.setDefaultProperties();
matRad_cfg.logLevel = 3;
cd(back);
clear back;
matRad_cfg.dispInfo('Restored default properties and returned to original directory!\n');