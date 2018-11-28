function [resultGUI,info] = matRad_directApertureOptimization(dij,cst,apertureInfo,resultGUI,pln,stf)
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
    addpath(fullfile(matRadRootDir,'tools'))
    [env, ~] = matRad_getEnvironment();
    
    switch env
         case 'MATLAB'
            % get handle to Matlab command window
            mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
            cw          = mde.getClient('Command Window');
            xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
            h_cw        = handle(xCmdWndView,'CallbackProperties');

            % set Key Pressed Callback of Matlab command window
            set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);
    end

end

% Set the IPOPT options.
matRad_ipoptOptions;

% set optimization options
options.radMod          = pln.radiationMode;
options.bioOpt          = pln.propOpt.bioOptimization;
options.ID              = [pln.radiationMode '_' pln.propOpt.bioOptimization];
options.FMO             = false; % let optimizer know that this is FMO
options.numOfScenarios  = dij.numOfScenarios;

% initialize global variables for optimizer
global matRad_global_x;
global matRad_global_d;
global matRad_STRG_C_Pressed;
global matRad_objective_function_value;
global matRad_global_apertureInfo;

matRad_global_x                 = NaN * ones(dij.totalNumOfBixels,1); % works with bixel weights even though we do dao!
matRad_global_d                 = NaN * ones(dij.numOfVoxels,1);
matRad_STRG_C_Pressed           = false;
matRad_objective_function_value = [];

% adjust overlap priorities
cst_Over = matRad_setOverlapPriorities(cst);

% adjust objectives _and_ constraints internally for fractionation 
for i = 1:size(cst_Over,1)
    for j = 1:size(cst_Over{i,6},1)
       cst_Over{i,6}(j).dose = cst_Over{i,6}(j).dose/pln.numOfFractions;
    end
end

if isfield(apertureInfo,'scaleFacRx')
    %weights were scaled to acheive 95% PTV coverage
    %scale back to "optimal" weights
    apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) = apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes)/apertureInfo.scaleFacRx;
end

if pln.propOpt.preconditioner
    %rescale dij matrix, so that apertureWeight/bixelWidth ~= 1
    % gradient wrt weights ~ 1, gradient wrt leaf pos
    % ~ apertureWeight/(bixelWidth) ~1
    
    % need to get the actual weights, so use the jacobiScale vector to
    % convert from the variables
    dij.scaleFactor = mean(apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes)./apertureInfo.jacobiScale)/(apertureInfo.bixelWidth);
    
    dij.weightToMU = dij.weightToMU*dij.scaleFactor;
    apertureInfo.weightToMU = apertureInfo.weightToMU*dij.scaleFactor;
    apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) = apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes)/dij.scaleFactor;
end

% update aperture info vector
apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);
apertureInfo.newIteration = true;
% define apertureInfo as a global vector to be updated once each iteration
matRad_global_apertureInfo = apertureInfo;

% set bounds on optimization variables
options.lb              = apertureInfo.limMx(:,1);                                          % Lower bound on the variables.
options.ub              = apertureInfo.limMx(:,2);                                          % Upper bound on the variables.
options.runVMAT         = pln.propOpt.runVMAT;
[options.cl,options.cu] = matRad_daoGetConstBounds(cst_Over,apertureInfo,options);   % Lower and upper bounds on the constraint functions.

% set callback functions.
funcs.objective         = @(x) matRad_daoObjFunc(x,dij,cst_Over,options);
funcs.iterfunc          = @(iter,objective,parameter) matRad_IpoptIterFunc(iter,objective,parameter,options.ipopt.max_iter);
funcs.gradient          = @(x) matRad_daoGradFunc(x,dij,cst_Over,options);
funcs.constraints       = @(x) matRad_daoConstFunc(x,dij,cst_Over,options);
funcs.jacobian          = @(x) matRad_daoJacobFunc(x,dij,cst_Over,options);
funcs.jacobianstructure = @( ) matRad_daoGetJacobStruct(apertureInfo,dij,cst_Over,options);

% Run IPOPT.
[optApertureInfoVec, info] = ipopt(apertureInfo.apertureVector,funcs,options);

% unset Key Pressed Callback of Matlab command window and delete waitbar
if ~isdeployed && strcmp(env,'MATLAB')
    set(h_cw, 'KeyPressedCallback',' ');
end

% clear global variables after optimization
switch env
    case 'MATLAB'
        clearvars -global matRad_global_x matRad_global_d;
    case 'OCTAVE' 
        clear -global matRad_global_x matRad_global_d;
end

if pln.propOpt.preconditioner
    % revert scaling
    
    dij.weightToMU = dij.weightToMU./dij.scaleFactor;
    resultGUI.apertureInfo.weightToMU = resultGUI.apertureInfo.weightToMU./dij.scaleFactor;
    optApertureInfoVec(1:apertureInfo.totalNumOfShapes) = optApertureInfoVec(1:apertureInfo.totalNumOfShapes).*dij.scaleFactor;
end

% update the apertureInfoStruct and calculate bixel weights
resultGUI.apertureInfo = matRad_daoVec2ApertureInfo(resultGUI.apertureInfo,optApertureInfoVec);

% override also bixel weight vector in optResult struct
resultGUI.w    = resultGUI.apertureInfo.bixelWeights;
resultGUI.wDao = resultGUI.apertureInfo.bixelWeights;

dij.scaleFactor = 1;

resultGUI.apertureInfo = matRad_preconditionFactors(resultGUI.apertureInfo);

% calc dose and reshape from 1D vector to 3D array
d = matRad_backProjection(resultGUI.w,dij,options);
resultGUI.physicalDose = reshape(d{1},dij.dimensions);

if isfield(pln,'scaleDRx') && pln.scaleDRx
    %Scale D95 in target to RXDose
    resultGUI.QI = matRad_calcQualityIndicators(cst,pln,resultGUI.physicalDose);
    
    resultGUI.apertureInfo.scaleFacRx = max((pln.DRx/pln.numOfFractions)./[resultGUI.QI(pln.RxStruct).D_95]');
    resultGUI.apertureInfo.apertureVector(1:resultGUI.apertureInfo.totalNumOfShapes) = resultGUI.apertureInfo.apertureVector(1:resultGUI.apertureInfo.totalNumOfShapes)*resultGUI.apertureInfo.scaleFacRx;
    
    % update the apertureInfoStruct and calculate bixel weights
    resultGUI.apertureInfo = matRad_daoVec2ApertureInfo(resultGUI.apertureInfo,resultGUI.apertureInfo.apertureVector);
    
    % override also bixel weight vector in optResult struct
    resultGUI.w    = resultGUI.apertureInfo.bixelWeights;
    resultGUI.wDao = resultGUI.apertureInfo.bixelWeights;
    
    resultGUI.physicalDose = resultGUI.physicalDose.*resultGUI.apertureInfo.scaleFacRx;
end

% update apertureInfoStruct with the maximum leaf speeds per segment
if pln.propOpt.runVMAT
    resultGUI.apertureInfo = matRad_maxLeafSpeed(resultGUI.apertureInfo);
    
    %optimize delivery
    resultGUI = matRad_optDelivery(resultGUI,1);
    resultGUI = matRad_calcDeliveryMetrics(resultGUI,pln,stf);
end

