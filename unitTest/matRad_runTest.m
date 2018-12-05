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
% limiting the cutoffLevel and lateralCutoff for faster computation
matRad_unitTestTextManipulation('matRad_calcPhotonDose.m', 'lateralCutoff = 50', 'lateralCutoff = 20;')
matRad_unitTestTextManipulation('matRad_calcParticleDose.m', 'cutOffLevel          = 0.99', '       cutOffLevel          = 0.8;')
matRad_unitTestTextManipulation('matRad_ipoptOptions.m', 'options.ipopt.max_iter', 'options.ipopt.max_iter = 10;', '../optimization/')

% supressing the inherent Ocatave warnings for division by zero
if strcmp(matRad_getEnvironment,'OCTAVE')
    warning("off", "Octave:divide-by-zero")
end

unitTestBool = true;

disp('Unit test run example 1');
matRad_example1_phantom
disp('Unit test run example 2');
matRad_example2_photons
disp('Unit test run example 3');
matRad_example3_photonsDAO
disp('Unit test run example 5');
matRad_example5_protons
disp('Unit test run example 6');
matRad_example6_protonsNoise
disp('Unit test run example 7');
matRad_example7_carbon
disp('Unit test run example 8');
matRad_example8_protonsRobust
disp('Unit test run matRad script');
matRad