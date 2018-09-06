% function statusAll = runMatRadTests(varargin)
%% This file runs the complete matRad test suite.

%% Set path
clc; clear all; close all
addpath(fullfile(pwd,'..'));
addpath(fullfile(pwd,'..','dicomImport'));
addpath(fullfile(pwd,'..','tools'));




status = []
status = test("matRad_example1_phantom")`
% status = [status; test("matRad_example2_photons")]
% status = [status; test("matRad_example5_protons")]

if any(~status)
  error ('test failure')
  else
  disp('All test suits has been passed!')
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