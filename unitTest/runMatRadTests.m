function statusAll = runMatRadTests(varargin)
%% This file runs the complete matRad test suite.

%% Set path

addpath(fullfile(pwd,'..'));
addpath(fullfile(pwd,'..','dicomImport'));
addpath(fullfile(pwd,'..','tools'));

matRad_example5_protons
matRad_example1_phantom
matRad_example2_photons

exit