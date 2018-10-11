%% This file runs the complete matRad test suite.

%% Set path
clc; clear; close all
addpath(genpath(fullfile(pwd,'..')));


% limiting the optimization to 10 iterations
cd ./../optimization/
fileId = fopen('matRad_ipoptOptions.m', 'a');
fprintf(fileId, '\n');
fwrite(fileId, 'options.ipopt.max_iter = 10;');
fclose(fileId);
cd ./../unitTest/

% supressing the inherent Ocatave warnings for division by zero
warning("off", "Octave:divide-by-zero")

param.logLevel = 3;

% matRad_example1_phantom
% matRad_example2_photons
% matRad_example3_photonsDAO
% matRad_example5_protons
% matRad_example6_protonsNoise
% matRad_example7_carbon
% matRad
matRad_example8_protonsRobust