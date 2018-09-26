%% This file runs the complete matRad test suite.

%% Set path
clc; clear all; close all
addpath(fullfile(pwd,'..'));

unitTests = {@matRad_example1_phantom, ...
             @matRad_example2_photons, ...
             @matRad_example3_photonsDAO, ...
             @matRad_example5_protons, ...
             @matRad_example6_protonsNoise, ...
             @matRad_example7_carbon, ...
             @matRad_example8_protonsRobust};

status = [];

for k = 1:length(unitTests);
    
    [cst, stf, pln, ct, dij, resultGUI] = unitTests{k}();

end