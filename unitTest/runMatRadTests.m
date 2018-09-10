% function statusAll = runMatRadTests(varargin)
%% This file runs the complete matRad test suite.

%% Set path
clc; clear all; close all
addpath(fullfile(pwd,'..'));
addpath(fullfile(pwd,'..','dicomImport'));
addpath(fullfile(pwd,'..','tools'));

[cst, stf, pln, ct] = matRad_example1_phantom();
report = matRad_unitTest_sizeCheck(cst);

status = [];
status = [status; report.status];

% status = []
% status = test("matRad_example1_phantom")
% % status = [status; test("matRad_example2_photons")]
% % status = [status; test("matRad_example5_protons")]

if any(~status)
  error ('Test suit has not passed, look into the warnings to find the source of the error.')
  else
  disp('All tests has been passed!')
end


% matRad_example1_phantom
% disp('example 1 done')
% matRad_example2_photons
% disp('example 2 done')
% matRad_example3_photonsDAO
% disp('example 3 done')
% matRad_example5_protons
% disp('example 5 done')
% matRad_example6_protonsNoise
% disp('example 6 done')
% matRad_example7_carbon
% disp('example 7 done')
% matRad_example8_protonsRobust
% disp('example 8 done')

% exit