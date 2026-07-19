function [resultGUI, optimizer] = matRad_directApertureOptimization(dij, cst, apertureInfo, resultGUI, pln)
% matRad function to run direct aperture optimization
%
% call:
%   [optResult,optimizer] = matRad_directApertureOptimization(dij,cst,apertureInfo,optResult,pln)
%
% input:
%   dij:            matRad dij struct
%   cst:            matRad cst struct
%   apertureInfo:   aperture shape info struct
%   optResult:      resultGUI struct to which the output data will be added, if
%                   this field is empty optResult struct will be created
%
%   pln:            matRad pln struct
%
% output:
%   optResult:  struct containing optimized fluence vector, dose, and
%               shape info
%   optimizer:  used optimizer object
%
% References
%   [1] http://dx.doi.org/10.1118/1.4914863
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

% promote an empty/incomplete resultGUI and seed its apertureInfo from the
% (always present) apertureInfo argument, so callers don't have to prime
% resultGUI.apertureInfo themselves before the first DAO call
if ~isfield(resultGUI, 'apertureInfo')
    resultGUI.apertureInfo = apertureInfo;
end

% adjust overlap priorities
cst = matRad_setOverlapPriorities(cst);

% check & adjust objectives and constraints internally for fractionation
for i = 1:size(cst, 1)
    for j = 1:numel(cst{i, 6})
        obj = cst{i, 6}{j};

        % In case it is a default saved struct, convert to object
        % Also intrinsically checks that we have a valid optimization
        % objective or constraint function in the end
        if ~isa(obj, 'matRad_DoseOptimizationFunction')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!', i, j);
            end
        end

        obj = obj.setDoseParameters(obj.getDoseParameters() / pln.numOfFractions);

        cst{i, 6}{j} = obj;
    end
end

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst, dij.ctGrid.x, dij.ctGrid.y, dij.ctGrid.z, ...
                             dij.doseGrid.x, dij.doseGrid.y, dij.doseGrid.z);

if ~isfield(pln, 'bioModel')
    pln.bioModel = 'none';
end

if ~isa(pln.bioModel, 'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel, pln.radiationMode);
end

% set optimization options
options.ixForOpt     = 1;
options.numOfScen    = 1;
options.scenProb     = 1;
options.quantityOpt  = pln.propOpt.quantityOpt;
options.model        = pln.bioModel.model;

% update aperture info vector
if isfield(apertureInfo, 'scaleFacRx')
    % weights were scaled to achieve 95% PTV coverage
    % scale back to "optimal" weights
    apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) = ...
        apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) / apertureInfo.scaleFacRx;
end

if ~isfield(pln.propOpt, 'runVMAT')
    pln.propOpt.runVMAT = false;
end

% the preconditioner flag travels with the apertureInfo (set by the
% sequencer from pln.propSeq.preconditioner)
preconditioner = matRad_getFieldOrDefault(apertureInfo, 'preconditioner', false);

if preconditioner
    % rescale dij matrix, so that apertureWeight/bixelWidth ~= 1
    % gradient wrt weights ~ 1, gradient wrt leaf pos
    % ~ apertureWeight/(bixelWidth) ~1

    % need to get the actual weights, so use the jacobiScale vector to
    % convert from the variables
    dij.scaleFactor = mean(apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) ./ apertureInfo.jacobiScale) / (apertureInfo.bixelWidth);

    dij.weightToMU = dij.weightToMU * dij.scaleFactor;
    apertureInfo.weightToMU = apertureInfo.weightToMU * dij.scaleFactor;
    apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) = apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) / dij.scaleFactor;
end

% Use Dose Projection only
backProjection = matRad_DoseProjection();

% the VMAT/DAO mode choice is made exactly once here; everything below
% dispatches polymorphically through the optiProb instance
if pln.propOpt.runVMAT
    optiProb = matRad_OptimizationProblemVMAT(backProjection, apertureInfo);
else
    optiProb = matRad_OptimizationProblemDAO(backProjection, apertureInfo);
end
apertureInfo = optiProb.matRad_daoVec2ApertureInfo(apertureInfo, apertureInfo.apertureVector);
apertureInfo.newIteration = true; % do we need this?
optiProb.apertureInfo = apertureInfo;

if ~isfield(pln.propOpt, 'optimizer')
    pln.propOpt.optimizer = 'IPOPT';
end

switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    otherwise
        matRad_cfg.dispWarning('Optimizer ''%s'' not known! Fallback to IPOPT!', pln.propOpt.optimizer);
        optimizer = matRad_OptimizerIPOPT;
end

% Run IPOPT.
optimizer = optimizer.optimize(apertureInfo.apertureVector, optiProb, dij, cst);
optApertureInfoVec = optimizer.wResult;

% Additional VMAT stuff
if preconditioner
    % revert the dij/vector scaling; resultGUI.apertureInfo (the template
    % for the conversion below) was never scaled, so its weightToMU
    % already holds the original value
    dij.weightToMU = dij.weightToMU ./ dij.scaleFactor;
    optApertureInfoVec(1:apertureInfo.totalNumOfShapes) = optApertureInfoVec(1:apertureInfo.totalNumOfShapes) .* dij.scaleFactor;
end

% update the apertureInfoStruct and calculate bixel weights
% Use optiProb dispatch to automatically choose VMAT / DAO code
newApertureInfo = optiProb.matRad_daoVec2ApertureInfo(resultGUI.apertureInfo, optApertureInfoVec);

% override also bixel weight vector in optResult struct
w    = newApertureInfo.bixelWeights;
wDao = newApertureInfo.bixelWeights;

dij.scaleFactor = 1;

newApertureInfo = matRad_preconditionFactors(newApertureInfo);

% logging final results
matRad_cfg.dispInfo('Calculating final cubes...\n');

% merge the computed dose cubes into the passed-in resultGUI struct (an
% empty input is promoted to a fresh struct) instead of overwriting it, so
% that pre-existing fields such as resultGUI.sequencing are preserved
doseCubes = matRad_calcCubes(w, dij);
fNames = fieldnames(doseCubes);
for f = 1:numel(fNames)
    resultGUI.(fNames{f}) = doseCubes.(fNames{f});
end
resultGUI.w    = w;
resultGUI.wDAO = wDao;
resultGUI.apertureInfo = newApertureInfo;
if isfield(resultGUI, 'sequencing') && isstruct(resultGUI.sequencing)
    resultGUI.sequencing.apertureInfo = newApertureInfo;
end

if isfield(pln, 'scaleDRx') && pln.scaleDRx
    % Scale D95 in target to RXDose
    resultGUI.QI = matRad_calcQualityIndicators(cst, pln, resultGUI.physicalDose);

    resultGUI.apertureInfo.scaleFacRx = max((pln.DRx / pln.numOfFractions) ./ [resultGUI.QI(pln.RxStruct).D_95]');
    resultGUI.apertureInfo.apertureVector(1:resultGUI.apertureInfo.totalNumOfShapes) = ...
        resultGUI.apertureInfo.apertureVector(1:resultGUI.apertureInfo.totalNumOfShapes) * resultGUI.apertureInfo.scaleFacRx;

    % update the apertureInfoStruct and calculate bixel weights
    resultGUI.apertureInfo = optiProb.matRad_daoVec2ApertureInfo(resultGUI.apertureInfo, resultGUI.apertureInfo.apertureVector);

    % override also bixel weight vector in optResult struct
    resultGUI.w    = resultGUI.apertureInfo.bixelWeights;
    resultGUI.wDao = resultGUI.apertureInfo.bixelWeights;

    resultGUI.physicalDose = resultGUI.physicalDose .* resultGUI.apertureInfo.scaleFacRx;
end

% update apertureInfoStruct with the maximum leaf speeds per segment
if pln.propOpt.runVMAT
    resultGUI.apertureInfo = matRad_OptimizationProblemVMAT.maxLeafSpeed(resultGUI.apertureInfo);

    % optimize delivery
    resultGUI.apertureInfo = matRad_OptimizationProblemVMAT.optDelivery(resultGUI.apertureInfo, 1);
end
