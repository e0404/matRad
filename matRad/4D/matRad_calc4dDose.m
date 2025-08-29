function resultGUI = matRad_calc4dDose(dij, pln,stf,resultGUI, totalPhaseMatrix)
% wrapper for the whole 4D dose calculation pipeline and calculated dose
% accumulation
%
% call
%   ct = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI)
%
% input
%   dij:             matRad dij struct
%   resultGUI:       struct containing optimized fluence vector
%   totalPhaseMatrix optional intput for totalPhaseMatrix
% output
%   resultGUI:      structure containing phase dose, RBE weighted dose, etc
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2025 the matRad development team.
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

if  ~exist('totalPhaseMatrix','var')
    sequencer = matRad_SequencingBase.getSequencerFromPln(pln);
    sequence = sequencer.sequence(resultGUI.w,stf);
    sequence = sequencer.makePhaseMatrix(sequence,pln.multScen.numOfCtScen, pln.multScen.motionPeriod);
    resultGUI.sequencing   = sequence;
    totalPhaseMatrix  = vertcat(sequence.phaseMatrix);
end


if ~isa(pln.bioModel,'matRad_BiologicalModel')
    pln.bioModel = matRad_BiologicalModel.validate(pln.bioModel,pln.radiationMode);
end

% compute all phases
for i = 1:size(totalPhaseMatrix,2)

    tmpResultGUI = matRad_calcCubes(totalPhaseMatrix(:,i),dij,i);

    % compute physical dose for physical opt
    resultGUI.phaseDose{i} = tmpResultGUI.physicalDose;
    if ~isa(pln.bioModel,'matRad_EmptyBiologicalModel')
        % compute RBExDose
        if isa(pln.bioModel,'matRad_ConstantRBE')
            resultGUI.phaseRBExDose{i} = tmpResultGUI.RBExDose;
        elseif isa(pln.bioModel,'matRad_LQBasedModel')
            resultGUI.phaseAlphaDose{i}    = tmpResultGUI.alpha .* tmpResultGUI.physicalDose;
            resultGUI.phaseSqrtBetaDose{i} = sqrt(tmpResultGUI.beta) .* tmpResultGUI.physicalDose;
            ix = dij.ax{i} ~=0;
            resultGUI.phaseEffect{i}    = resultGUI.phaseAlphaDose{i} + resultGUI.phaseSqrtBetaDose{i}.^2;
            resultGUI.phaseRBExDose{i}     = zeros(ct.cubeDim);
            resultGUI.phaseRBExDose{i}(ix) = ((sqrt(dij.ax{i,1}(ix).^2 + 4 .* dij.bx{i,1}(ix) .* resultGUI.phaseEffect{i}(ix)) -dij.ax{i,1}(ix))./(2.*dij.bx{i,1}(ix)));
        else
            matRad_cfg.dispError('Unsupported biological model %s!',pln.bioModel.model);
        end
    end

     for beamIx = 1:dij.numOfBeams
        resultGUI.(['phaseDose_beam', num2str(beamIx)]){i} = tmpResultGUI.(['physicalDose_beam', num2str(beamIx)]);
        if ~isa(pln.bioModel,'matRad_EmptyBiologicalModel')
            % compute RBExD 
            if isa(pln.bioModel,'matRad_ConstantRBE')
                resultGUI.(['phaseRBExDose_beam', num2str(beamIx)]){i} = tmpResultGUI.(['RBExDose_beam', num2str(beamIx)]);
            elseif isa(pln.bioModel,'matRad_LQBasedModel')
                resultGUI.(['phaseAlphaDose_beam', num2str(beamIx)]){i} = tmpResultGUI.(['alpha_beam', num2str(beamIx)]).*tmpResultGUI.(['physicalDose_beam', num2str(beamIx)]);
                resultGUI.(['phaseSqrtBetaDose_beam', num2str(beamIx)]){i} = sqrt(tmpResultGUI.(['beta_beam', num2str(beamIx)])).*tmpResultGUI.(['physicalDose_beam', num2str(beamIx)]);
                 ix = dij.ax{i} ~=0;
                 resultGUI.(['phaseEffect_beam', num2str(beamIx)]){i}    = resultGUI.(['phaseAlphaDose_beam', num2str(beamIx)]){i} + resultGUI.(['phaseSqrtBetaDose_beam', num2str(beamIx)]){i}.^2;  
                 resultGUI.(['phaseRBExDose_beam', num2str(beamIx)]){i}    = zeros(ct.cubeDim);
                 resultGUI.(['phaseRBExDose_beam', num2str(beamIx)]){i}     =  ((sqrt(dij.ax{i,1}(ix).^2 + 4 .* dij.bx{i,1}(ix) .* resultGUI.(['phaseEffect_beam', num2str(beamIx)]){i}(ix)) - dij.ax{i,1}(ix))./(2.*dij.bx{i,1}(ix)));
            else
                matRad_cfg.dispError('Unsupported biological model %s!',pln.bioModel.model);
            end
        end
    end
end



end

