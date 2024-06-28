function result = matRad_runTests(folder,withCoverage)
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

back = cd(fileparts(mfilename("fullpath")));
matRad_cfg = matRad_rc;

if nargin < 2
    withCoverage = false;
end

if withCoverage && matRad_cfg.isOctave
    matRad_cfg.dispWarning('Coverage collection not possible with Octave. Turning off!');
    withCoverage = false;
end

%add MOxUnit submodule if not yet on path
if exist('moxunit_runtests','file') == 2
    matRad_cfg.dispInfo('MOxUnit is already on the path and this version will be used for testing: ');
else
    matRad_cfg.dispInfo('MOxUnit submodule will be added to the path and used for testing: ');
    rootFolder = cd(fullfile(matRad_cfg.matRadRoot,'submodules','MOxUnit','MOxUnit'));
    moxunit_set_path;
    cd(rootFolder);
end
matRad_cfg.dispInfo('%s\n',fileparts(which("moxunit_runtests")));

%add MOcov if coverage will be recorded
if withCoverage
    if exist('mocov','file') == 2
        matRad_cfg.dispInfo('MOcov is already on the path and this version will be used for coverage collection: ');
    else
        matRad_cfg.dispInfo('MOcov submodule will be added to the path and used for coverage collection: ');
        addpath(fullfile(matRad_cfg.matRadRoot,'submodules','MOcov','MOcov'));
    end
    matRad_cfg.dispInfo('%s\n',fileparts(which("mocov.m")));
end


matRad_cfg.dispInfo('Setting default properties for testing and starting test suite!\n');
matRad_cfg.setDefaultPropertiesForTesting();
matRad_cfg.logLevel = 1;

addpath(fullfile(matRad_cfg.matRadRoot,'test'));
if nargin < 1
    folder = 'test';
end

%add test folder to path
addpath(genpath(fullfile(matRad_cfg.matRadRoot,folder)));

if withCoverage
    result = moxunit_runtests(folder,'-recursive','-junit_xml_file','testresults.xml',...
        '-with_coverage','-cover','matRad',...
        '-cover_xml_file','coverage.xml','-cover_json_file','coverage.json',...
        '-cover_method','profile');
else
    result = moxunit_runtests(folder,'-recursive','-junit_xml_file','testresults.xml');
end


matRad_cfg.setDefaultProperties();
matRad_cfg.logLevel = 3;
cd(back);
clear back;
matRad_cfg.dispInfo('Restored default properties and returned to original directory!\n');
