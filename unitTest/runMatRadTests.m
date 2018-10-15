%% This file runs the complete matRad test suite.

%% Set path
clc; clear; close all
addpath(genpath(fullfile(pwd,'..')));


% limiting the optimization to 10 iterations for faster computation
% limiting the cutoffLevel and lateralCutoff for faster computation
matRad_unitTestTextManipulation()

% supressing the inherent Ocatave warnings for division by zero
warning("off", "Octave:divide-by-zero")

param.logLevel = 3;

matRad_example1_phantom
matRad_example2_photons
matRad_example3_photonsDAO
matRad_example5_protons
matRad_example6_protonsNoise
matRad_example7_carbon
matRad_example8_protonsRobust
matRad