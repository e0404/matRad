%% This file runs the complete matRad test suite.
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

%% Set path
run(['..' filesep 'matRad_rc']);

%% Prepare settings for testing
global matRad_cfg; matRad_cfg = MatRad_Config.instance();
matRad_cfg.setDefaultPropertiesForTesting();

% supressing the inherent Ocatave warnings for division by zero
if strcmp(matRad_getEnvironment,'OCTAVE')
    warning('off','Octave:divide-by-zero');
end

exampleScripts = {'matRad_example1_phantom.m',...
    'matRad_example2_photons.m',...
    'matRad_example3_photonsDAO.m',...
    'matRad_example4_photonsMC.m',...
    'matRad_example5_protons.m',...
    'matRad_example6_protonsNoise.m',...
    'matRad_example7_carbon.m'};

unitTestBixelWidth = 20;

matRad_unitTestTextManipulation(exampleScripts,'pln.propStf.bixelWidth',['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth)], '../examples/');
matRad_unitTestTextManipulation(exampleScripts,'display(','%%%%%%%%%%%%%%% REMOVED DISPLAY FOR TESTING %%%%%%%%%%%%%%', '../examples/');
matRad_unitTestTextManipulation('matRad.m','pln.propStf.bixelWidth',['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth)], '../');

disp('Integration test run example 1');
matRad_example1_phantom
disp('Integration test run example 2');
matRad_example2_photons
disp('Integration test run example 3');
matRad_example3_photonsDAO
disp('Integration test run example 4');
matRad_example4_photonsMC
disp('Integration test run example 5');
matRad_example5_protons
disp('Integration test run example 6');
matRad_example6_protonsNoise
disp('Integration test run example 7');
matRad_example7_carbon
disp('Integration test run matRad script');
matRad