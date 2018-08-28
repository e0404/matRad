function statusAll = runMatRadTests(varargin)
%% This file runs the complete matRad test suite.

%% Set path

addpath(fullfile(pwd,'..'));
addpath(fullfile(pwd,'..','dicomImport'));
addpath(fullfile(pwd,'..','tools'));

matRad_example5_protons
disp('example 5 done')
matRad_example1_phantom
disp('example 1 done')
matRad_example2_photons
disp('example 2 done')
matRad_example3_photonsDAO
disp('example 3 done')

exit