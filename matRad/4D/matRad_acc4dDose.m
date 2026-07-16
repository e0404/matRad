function resultGUI = matRad_acc4dDose( dij, pln, ct, cst,resultGUI, accType)
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
matRad_cfg = MatRad_Config.instance();

if ~isa(pln.bioModel,'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel,pln.radiationMode);
end

% accumulation
resultGUI.accPhysicalDose = matRad_doseAcc(ct,resultGUI.phaseDose, cst, accType);
if isa(pln.bioModel,'matRad_ConstantRBE')

    resultGUI.accRBExDose = matRad_doseAcc(ct,resultGUI.phaseRBExDose, cst, accType);

elseif isa(pln.bioModel,'matRad_LQBasedModel')

    resultGUI.accAlphaDose    = matRad_doseAcc(ct,resultGUI.phaseAlphaDose, cst,accType);
    resultGUI.accSqrtBetaDose = matRad_doseAcc(ct,resultGUI.phaseSqrtBetaDose, cst, accType);

    % only compute where we have biologically defined tissue
    ix = (ax{1} ~= 0);

    resultGUI.accEffect = resultGUI.accAlphaDose + resultGUI.accSqrtBetaDose.^2;

    resultGUI.accRBExDose     = zeros(ct.cubeDim);
    resultGUI.accRBExDose(ix) = ((sqrt(ax{1}(ix).^2 + 4 .* bx{1}(ix) .* resultGUI.accEffect(ix)) - ax{1}(ix))./(2.*bx{1}(ix)));
end

for beamIx = 1:dij.numOfBeams
    resultGUI.(['accPhysicalDose_beam', num2str(beamIx)])= matRad_doseAcc(ct,resultGUI.(['phaseDose_beam', num2str(beamIx)]), cst, accType);
    if isa(pln.bioModel,'matRad_ConstantRBE')
        resultGUI.(['accRBExDose_beam', num2str(beamIx)]) = matRad_doseAcc(ct,resultGUI.(['phaseRBExDose_beam', num2str(beamIx)]), cst, accType);
   elseif isa(pln.bioModel,'matRad_LQBasedModel')
        resultGUI.(['accAlphaDose_beam', num2str(beamIx)])      = matRad_doseAcc(ct,resultGUI.(['phaseAlphaDose_beam', num2str(beamIx)]), cst, accType);
        resultGUI.(['accSqrtBetaDose_beam', num2str(beamIx)])   = matRad_doseAcc(ct,resultGUI.(['phaseAlphaDose_beam', num2str(beamIx)]), cst, accType);
        resultGUI.(['accEffect_beam', num2str(beamIx)])         = resultGUI.(['accAlphaDose_beam', num2str(beamIx)]) + resultGUI.(['accSqrtBetaDose_beam', num2str(beamIx)]).^2;
        resultGUI.(['accRBExDose_beam', num2str(beamIx)]){i}       = zeros(ct.cubeDim);
        resultGUI.(['accRBExDose_beam', num2str(beamIx)]){i}       =  ((sqrt(ax{i}(ix).^2 + 4 .* bx{i}(ix) .* resultGUI.(['accEffect_beam', num2str(beamIx)]){i}(ix)) - ax{i}(ix))./(2.*bx{i}(ix)));    
    end
end


end