function matRad_runTests(folder,withCoverage)
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

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispInfo('Setting default properties for testing and starting test suite!\n');
matRad_cfg.setDefaultPropertiesForTesting();
matRad_cfg.logLevel = 1;

if nargin < 2
    withCoverage = false;
end

if withCoverage && matRad_cfg.isOctave
    matRad_cfg.dispWarning('Coverage collection not possible with Octave. Turning off!');
    withCoverage = false;
end

if nargin < 1
    folder = 'test';
end


back = cd(matRad_cfg.matRadRoot);

if withCoverage
    moxunit_runtests('test','-recursive','-junit_xml_file','testresults.xml',...
        '-with_coverage','-cover','.','-cover_xml_file','coverage.xml','-cover_json_file','coverage.json',...
        '-cover_exclude','submodules','-cover_exclude','examples','-cover_method','profile');
else
    moxunit_runtests(folder,'-recursive');
end


matRad_cfg.setDefaultProperties();
matRad_cfg.logLevel = 3;
cd(back);
clear back;
matRad_cfg.dispInfo('Restored default properties and returned to original directory!\n');
