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
run(['..' filesep 'matRad_rc'])


% limiting the optimization to 10 iterations for faster computation
matRad_unitTestTextManipulation('matRad_OptimizerIPOPT.m','obj.options.max_iter','obj.options.max_iter = 10;','../optimization/optimizer/');
% limiting the cutoffLevel and lateralCutoff for faster computation
matRad_unitTestTextManipulation('matRad_calcPhotonDose.m', 'lateralCutoff = ', 'lateralCutoff = 20;')
matRad_unitTestTextManipulation('matRad_calcParticleDose.m', 'cutOffLevel = ', '       cutOffLevel          = 0.8;')
% limit number of histories for MC to 100
matRad_unitTestTextManipulation('matRad_calcParticleDoseMC.m', '    nCasePerBixel', '    nCasePerBixel = 100;')
matRad_unitTestTextManipulation('matRad_calcPhotonDoseMC.m', '    nCasePerBixel', '    nCasePerBixel = 100;')
matRad_unitTestTextManipulation('matRad_calcDoseDirectMC.m', '  nHistories = 2e4;', '  nHistories = 100;')



exampleScripts = {'matRad_example1_phantom.m',...
    'matRad_example2_photons.m',...
    'matRad_example3_photonsDAO.m',...
    'matRad_example4_photonsMC.m',...
    'matRad_example5_protons.m',...
    'matRad_example6_protonsNoise.m',...
    'matRad_example7_carbon.m'};

unitTestBixelWidth = 20;

matRad_unitTestTextManipulation(exampleScripts,'pln.propStf.bixelWidth',['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth)], '../examples/');
matRad_unitTestTextManipulation(exampleScripts,'display(','%%%%%%%%%%%%%%% REMOVED DISPLAY FOR UNIT TESTING %%%%%%%%%%%%%%', '../examples/');
matRad_unitTestTextManipulation(exampleScripts,'matRadGUI','%%%%%%%%%%%%%%% REMOVED matRadGUI FOR UNIT TESTING %%%%%%%%%%%%%%', '../examples/');
matRad_unitTestTextManipulation('matRad.m','matRadGUI','%%%%%%%%%%%%%%% REMOVED matRadGUI FOR UNIT TESTING %%%%%%%%%%%%%%', '../');
matRad_unitTestTextManipulation('matRad.m','pln.propStf.bixelWidth',['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth)], '../');

% set coarser and anisotropic dose grid for unit testing
doseCalcResX = 5;
doseCalcResY = 6;
doseCalcResZ = 7;
matRad_unitTestTextManipulation(exampleScripts,'pln.propDoseCalc.doseGrid.resolution.x',['pln.propDoseCalc.doseGrid.resolution.x = ' num2str(doseCalcResX)], '../examples/');
matRad_unitTestTextManipulation(exampleScripts,'pln.propDoseCalc.doseGrid.resolution.y',['pln.propDoseCalc.doseGrid.resolution.y = ' num2str(doseCalcResY)], '../examples/');
matRad_unitTestTextManipulation(exampleScripts,'pln.propDoseCalc.doseGrid.resolution.z',['pln.propDoseCalc.doseGrid.resolution.z = ' num2str(doseCalcResZ)], '../examples/');

% supressing the inherent Ocatave warnings for division by zero
if strcmp(matRad_getEnvironment,'OCTAVE')
    warning('off','Octave:divide-by-zero')
end

unitTestBool = true;

disp('Unit test run example 1');
matRad_example1_phantom
disp('Unit test run example 2');
matRad_example2_photons
disp('Unit test run example 3');
matRad_example3_photonsDAO
disp('Unit test run example 4');
matRad_example4_photonsMC
disp('Unit test run example 5');
matRad_example5_protons
disp('Unit test run example 6');
matRad_example6_protonsNoise
disp('Unit test run example 7');
matRad_example7_carbon
disp('Unit test run matRad script');
matRad