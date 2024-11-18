function [resultGUI, timeSequence] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, totalPhaseMatrix,accType)
% wrapper for the whole 4D dose calculation pipeline and calculated dose
% accumulation
%
% call
%   ct = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI)
%
% input
%   ct :             ct cube
%   pln:             matRad plan meta information struct
%   dij:             matRad dij struct
%   stf:             matRad steering information struct
%   cst:             matRad cst struct
%   resultGUI:       struct containing optimized fluence vector
%   totalPhaseMatrix optional intput for totalPhaseMatrix
%   accType:         witch algorithim for dose accumulation
% output
%   resultGUI:      structure containing phase dose, RBE weighted dose, etc
%   timeSequence:   timing information about the irradiation
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('accType','var')
    accType = 'DDM';
end

if ~exist('totalPhaseMatrix','var')
    % make a time sequence for when each bixel is irradiated, the sequence
    % follows the backforth spot scanning
    timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI);

    % prepare a phase matrix
    motion       = 'linear'; % the assumed motion type
    timeSequence = matRad_makePhaseMatrix(timeSequence, ct.numOfCtScen, ct.motionPeriod, motion);

    resultGUI.bioModel = pln.bioModel;

    % the total phase matrix determines what beamlet will be administered in what ct phase
    totalPhaseMatrix   = vertcat(timeSequence.phaseMatrix);
else
    timeSequence = [];

end

% Get Biological Model
if ~isfield(pln,'bioModel')
    pln.bioModel = 'none';
end
if ~isa(pln.bioModel,'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel,pln.radiationMode);
end

if isa(pln.bioModel,'matRad_LQBasedModel')
    [ax,bx] = matRad_getPhotonLQMParameters(cst,numel(resultGUI.physicalDose));
end

% compute all phases
for i = 1:ct.numOfCtScen

    tmpResultGUI = matRad_calcCubes(totalPhaseMatrix(:,i),dij,i);

    % compute physical dose for physical opt
    if isa(pln.bioModel,'matRad_EmptyBiologicalModel')
        resultGUI.phaseDose{i} = tmpResultGUI.physicalDose;
        % compute RBExDose with const RBE
    elseif isa(pln.bioModel,'matRad_ConstantRBE')
        resultGUI.phaseRBExDose{i} = tmpResultGUI.RBExDose;
        % compute all fields
    elseif isa(pln.bioModel,'matRad_LQBasedModel')
        resultGUI.phaseAlphaDose{i}    = tmpResultGUI.alpha .* tmpResultGUI.physicalDose;
        resultGUI.phaseSqrtBetaDose{i} = sqrt(tmpResultGUI.beta) .* tmpResultGUI.physicalDose;
        ix = ax{i} ~=0;
        resultGUI.phaseEffect{i}    = resultGUI.phaseAlphaDose{i} + resultGUI.phaseSqrtBetaDose{i}.^2;
        resultGUI.phaseRBExDose{i}     = zeros(ct.cubeDim);
        resultGUI.phaseRBExDose{i}(ix) = ((sqrt(ax{i}(ix).^2 + 4 .* bx{i}(ix) .* resultGUI.phaseEffect{i}(ix)) - ax{i}(ix))./(2.*bx{i}(ix)));
    else
        matRad_cfg.dispError('Unsupported biological model %s!',pln.bioModel.model);
    end
end

% accumulation
if isa(pln.bioModel,'matRad_EmptyBiologicalModel')

    resultGUI.accPhysicalDose = matRad_doseAcc(ct,resultGUI.phaseDose, cst, accType);

elseif isa(pln.bioModel,'matRad_ConstantRBE')

    resultGUI.accRBExDose = matRad_doseAcc(ct,resultGUI.phaseRBExDose, cst, accType);

elseif isa(pln.bioModel,'matRad_LQBasedModel')

    resultGUI.accAlphaDose    = matRad_doseAcc(ct,resultGUI.phaseAlphaDose, cst,accType);
    resultGUI.accSqrtBetaDose = matRad_doseAcc(ct,resultGUI.phaseSqrtBetaDose, cst, accType);

    % only compute where we have biologically defined tissue
    ix = (ax{1} ~= 0);

    resultGUI.accEffect = resultGUI.accAlphaDose + resultGUI.accSqrtBetaDose.^2;

    resultGUI.accRBExDose     = zeros(ct.cubeDim);
    resultGUI.accRBExDose(ix) = ((sqrt(ax{1}(ix).^2 + 4 .* bx{1}(ix) .* resultGUI.accEffect(ix)) - ax{1}(ix))./(2.*bx{1}(ix)));
else
    matRad_cfg.dispError('Unsupported biological model %s!',pln.bioModel.model);
end
end

