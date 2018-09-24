% function statusAll = runMatRadTests(varargin)
%% This file runs the complete matRad test suite.

%% Set path
clc; clear all; close all
addpath(fullfile(pwd,'..'));

verified = load('verified.mat');

unitTests = {@matRad_example1_phantom, ...
             @matRad_example2_photons};

status = [];

for k = 1:length(unitTests);
    
    verVars = verified.(['example',num2str(k)]);
    
    [cst, stf, pln, ct, dij, resultGUI] = unitTests{k}();
    
    report = matRad_unitTest_sizeCheck(verVars, cst, ct, stf, pln, dij, resultGUI);
    
    status = [status; report.status];
    
    close all, clear cst stf pln ct dij resultGUI report
    
end


if any(~status)
  error ('Test suit didn''t pass!')
  else
  disp('All tests has been passed!')
end

