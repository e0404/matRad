%% This file runs the complete matRad test suite.

%% Set path
clc; clear; close all
addpath(genpath(fullfile(pwd,'..')));


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

matRad_example1_phantom
matRad_example2_photons
matRad_example3_photonsDAO
matRad_example5_protons
matRad_example6_protonsNoise
matRad_example7_carbon
matRad_example8_protonsRobust
matRad