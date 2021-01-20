function [optResult,optimizer] = matRad_directApertureOptimization(dij,cst,apertureInfo,optResult,pln)
% matRad function to run direct aperture optimization
%
% call
%   [optResult,optimizer] = matRad_directApertureOptimization(dij,cst,apertureInfo,pln)
%   [optResult,optimizer] = matRad_directApertureOptimization(dij,cst,apertureInfo,optResult,pln)
%
% input
%   dij:            matRad dij struct
%   cst:            matRad cst struct
%   apertureInfo:   aperture shape info struct
%   optResult:      resultGUI struct to which the output data will be added, if
%                   this field is empty optResult struct will be created
%                   (optional)
%   pln:            matRad pln struct
%
% output
%   optResult:  struct containing optimized fluence vector, dose, and
%               shape info
%   optimizer:  used optimizer object
%
% References
%   [1] http://dx.doi.org/10.1118/1.4914863
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


matRad_cfg = MatRad_Config.instance();

% adjust overlap priorities
cst = matRad_setOverlapPriorities(cst);

% check & adjust objectives and constraints internally for fractionation 
for i = 1:size(cst,1)
    for j = 1:numel(cst{i,6})
        obj = cst{i,6}{j};
        
        %In case it is a default saved struct, convert to object
        %Also intrinsically checks that we have a valid optimization
        %objective or constraint function in the end
        if ~isa(obj,'matRad_DoseOptimizationFunction')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
            end
        end
        
        obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
        
        cst{i,6}{j} = obj;        
    end
end

% resizing cst to dose cube resolution 
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                                 dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% update aperture info vector
apertureInfo = matRad_OptimizationProblemDAO.matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);

%Use Dose Projection only
backProjection = matRad_DoseProjection();

optiProb = matRad_OptimizationProblemDAO(backProjection,apertureInfo);

if ~isfield(pln.propOpt,'optimizer')
    pln.propOpt.optimizer = 'IPOPT';
end

switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    otherwise
        matRad_cfg.dispWarning('Optimizer ''%s'' not known! Fallback to IPOPT!',pln.propOpt.optimizer);
        optimizer = matRad_OptimizerIPOPT;
end

% Run IPOPT.
optimizer = optimizer.optimize(apertureInfo.apertureVector,optiProb,dij,cst);
wOpt = optimizer.wResult;

% update the apertureInfoStruct and calculate bixel weights
apertureInfo = matRad_OptimizationProblemDAO.matRad_daoVec2ApertureInfo(apertureInfo,wOpt);

% logging final results
matRad_cfg.dispInfo('Calculating final cubes...\n');
resultGUI = matRad_calcCubes(apertureInfo.bixelWeights,dij);
resultGUI.w    = apertureInfo.bixelWeights;
resultGUI.wDAO = apertureInfo.bixelWeights;
resultGUI.apertureInfo = apertureInfo;

