function [optResult,info] = matRad_directApertureOptimization(dij,cst,apertureInfo,optResult,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to run direct aperture optimization
%
% call
%   o[optResult,info] = matRad_directApertureOptimization(dij,cst,apertureInfo,optResult,pln,visBool)
%
% input
%   dij:            matRad dij struct
%   cst:            matRad cst struct
%   apertureInfo:   aperture shape info struct
%   optResult:      resultGUI struct to which the output data will be added, if
%                   this field is empty optResult struct will be created
%                   (optional)
%   pln:            matRad pln struct
%   visBool:        plots the objective function value in dependence of the
%                   number of iterations
%
% output
%   optResult:  struct containing optimized fluence vector, dose, and
%               shape info
%   info:       struct containing information about optimization
%
% References
%   [1] http://dx.doi.org/10.1118/1.4914863
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isdeployed % only if _not_ running as standalone
    
    % add path for optimization functions    
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'optimization'))
    
    % get handle to Matlab command window
    mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
    cw          = mde.getClient('Command Window');
    xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
    h_cw        = handle(xCmdWndView,'CallbackProperties');

    % set Key Pressed Callback of Matlab command window
    set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);

end

% initialize global variables for optimizer
global matRad_global_x;
global matRad_global_d;
global matRad_STRG_C_Pressed;
global matRad_objective_function_value;

matRad_global_x                 = NaN * ones(dij.totalNumOfBixels,1); % works with bixel weights even though we do dao!
matRad_global_d                 = NaN * ones(dij.numOfVoxels,1);
matRad_STRG_C_Pressed           = false;
matRad_objective_function_value = [];

% adjust overlap priorities
cst = matRad_setOverlapPriorities(cst);

% adjust objectives _and_ constraints internally for fractionation 
for i = 1:size(cst,1)
    for j = 1:size(cst{i,6},1)
       cst{i,6}(j).dose = cst{i,6}(j).dose/pln.numOfFractions;
    end
end

% update aperture info vector
apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);

% Set the IPOPT options.
matRad_ipoptOptions;

% set bounds on optimization variables
options.lb              = apertureInfo.limMx(:,1);                                          % Lower bound on the variables.
options.ub              = apertureInfo.limMx(:,2);                                          % Upper bound on the variables.
if isfield(pln,'VMAT') && pln.VMAT
    [options.cl,options.cu] = matRad_daoGetConstBounds(cst,apertureInfo,dij.numOfScenarios,pln.bioOptimization,pln.leafSpeedCst,pln.doseRateCst);   % Lower and upper bounds on the constraint functions.
else
    [options.cl,options.cu] = matRad_daoGetConstBounds(cst,apertureInfo,dij.numOfScenarios,pln.bioOptimization);   % Lower and upper bounds on the constraint functions.
end

% set optimization options
options.radMod          = pln.radiationMode;
options.bioOpt          = pln.bioOptimization;
options.ID              = [pln.radiationMode '_' pln.bioOptimization];
options.numOfScenarios  = dij.numOfScenarios;

% set callback functions.
funcs.objective         = @(x) matRad_daoObjFunc(x,apertureInfo,dij,cst,options);
funcs.gradient          = @(x) matRad_daoGradFunc(x,apertureInfo,dij,cst,options);
funcs.iterfunc          = @(iter,objective,paramter) matRad_IpoptIterFunc(iter,objective,paramter,options.ipopt.max_iter);
funcs.constraints       = @(x) matRad_daoConstFunc(x,apertureInfo,dij,cst,options);
funcs.jacobian          = @(x) matRad_daoJacobFunc(x,apertureInfo,dij,cst,options);
funcs.jacobianstructure = @( ) matRad_daoGetJacobStruct(apertureInfo,dij,cst);

% Run IPOPT.
[optApertureInfoVec, info] = ipopt(apertureInfo.apertureVector,funcs,options);

% unset Key Pressed Callback of Matlab command window and delete waitbar
if ~isdeployed
    set(h_cw, 'KeyPressedCallback',' ');
end

% clear global variables after optimization
clearvars -global matRad_global_x matRad_global_d;

% update the apertureInfoStruct and calculate bixel weights
optResult.apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,optApertureInfoVec);

% override also bixel weight vector in optResult struct
optResult.w    = optResult.apertureInfo.bixelWeights;
optResult.wDao = optResult.apertureInfo.bixelWeights;

% calc dose and reshape from 1D vector to 3D array
optResult.physicalDose = reshape(dij.physicalDose{1}*optResult.w,dij.dimensions);
