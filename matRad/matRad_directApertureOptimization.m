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
%   pln:            matRad pln struct. Relevant optional settings:
%                   pln.propOpt.runVMAT              optimize a VMAT arc instead
%                                                    of static apertures
%                   pln.propOpt.scaleToPrescription  scale the optimized plan
%                                                    such that the target D95
%                                                    reaches the prescription
%                   pln.propOpt.prescribedDose       total prescribed dose over
%                                                    all fractions [Gy]
%                   pln.propOpt.prescriptionStructIx cst row indices of the
%                                                    target structure(s) the
%                                                    prescription refers to (the
%                                                    worst D95 among them is
%                                                    scaled to the prescription)
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

% accept aperture info saved by an older matRad version
apertureInfo = matRad_upgradeApertureInfo(apertureInfo);

% promote an empty/incomplete resultGUI and seed its apertureInfo from the
% (always present) apertureInfo argument, so callers don't have to prime
% resultGUI.apertureInfo themselves before the first DAO call
if ~isfield(resultGUI, 'apertureInfo')
    resultGUI.apertureInfo = apertureInfo;
else
    resultGUI.apertureInfo = matRad_upgradeApertureInfo(resultGUI.apertureInfo);
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
if isfield(apertureInfo, 'prescriptionScaleFactor')
    % weights were scaled to reach the prescribed target coverage
    % scale back to "optimal" weights
    apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) = ...
        apertureInfo.apertureVector(1:apertureInfo.totalNumOfShapes) / apertureInfo.prescriptionScaleFactor;
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
resultGUI = matRad_mergeDoseCubes(resultGUI, matRad_calcCubes(w, dij));
resultGUI.w    = w;
resultGUI.wDAO = wDao;
resultGUI.apertureInfo = newApertureInfo;
if isfield(resultGUI, 'sequencing') && isstruct(resultGUI.sequencing)
    resultGUI.sequencing.apertureInfo = newApertureInfo;
end

% honor the pre-release top-level location of the prescription scaling
% settings for one release
if isfield(pln, 'scaleDRx')
    matRad_cfg.dispDeprecationWarning(['pln.scaleDRx/DRx/RxStruct are deprecated. Use pln.propOpt.scaleToPrescription/' ...
                                       'prescribedDose/prescriptionStructIx instead!']);
    pln.propOpt.scaleToPrescription  = pln.scaleDRx;
    pln.propOpt.prescribedDose       = matRad_getFieldOrDefault(pln, 'DRx', []);
    pln.propOpt.prescriptionStructIx = matRad_getFieldOrDefault(pln, 'RxStruct', []);
end

if matRad_getFieldOrDefault(pln.propOpt, 'scaleToPrescription', false)
    % scale the plan so that the worst target D95 matches the prescribed
    % dose per fraction
    if ~isfield(pln.propOpt, 'prescribedDose') || isempty(pln.propOpt.prescribedDose) || ...
       ~isfield(pln.propOpt, 'prescriptionStructIx') || isempty(pln.propOpt.prescriptionStructIx)
        matRad_cfg.dispError(['pln.propOpt.scaleToPrescription requires pln.propOpt.prescribedDose ' ...
                              'and pln.propOpt.prescriptionStructIx to be set!']);
    end

    % quality indicators of the unscaled plan (only used to derive the
    % scaling factor - the final QI are computed in matRad_planAnalysis)
    qi = matRad_calcQualityIndicators(cst, pln, resultGUI.physicalDose);

    prescriptionScaleFactor = max((pln.propOpt.prescribedDose / pln.numOfFractions) ./ ...
                                  [qi(pln.propOpt.prescriptionStructIx).D_95]');
    resultGUI.apertureInfo.prescriptionScaleFactor = prescriptionScaleFactor;
    resultGUI.apertureInfo.apertureVector(1:resultGUI.apertureInfo.totalNumOfShapes) = ...
        resultGUI.apertureInfo.apertureVector(1:resultGUI.apertureInfo.totalNumOfShapes) * prescriptionScaleFactor;

    % update the apertureInfo struct and recompute the bixel weights
    resultGUI.apertureInfo = optiProb.matRad_daoVec2ApertureInfo(resultGUI.apertureInfo, resultGUI.apertureInfo.apertureVector);

    resultGUI.w    = resultGUI.apertureInfo.bixelWeights;
    resultGUI.wDAO = resultGUI.apertureInfo.bixelWeights;

    % recompute all dose cubes from the scaled weights so that every
    % result quantity stays consistent with the scaled plan
    resultGUI = matRad_mergeDoseCubes(resultGUI, matRad_calcCubes(resultGUI.w, dij));
end

% update apertureInfoStruct with the maximum leaf speeds per segment
if pln.propOpt.runVMAT
    resultGUI.apertureInfo = matRad_OptimizationProblemVMAT.maxLeafSpeed(resultGUI.apertureInfo);

    % optimize delivery
    resultGUI.apertureInfo = matRad_OptimizationProblemVMAT.optDelivery(resultGUI.apertureInfo, 1);
end

end

function resultGUI = matRad_mergeDoseCubes(resultGUI, doseCubes)
% copy the dose cubes field-by-field into resultGUI, preserving any fields
% (e.g. resultGUI.sequencing) that matRad_calcCubes does not produce
fNames = fieldnames(doseCubes);
for f = 1:numel(fNames)
    resultGUI.(fNames{f}) = doseCubes.(fNames{f});
end
end
