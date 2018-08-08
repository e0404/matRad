function [resultGUI, timeSequence] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrapper for the whole 4D dose calculation pipeline and calculated dose
% aacumulation
%
% call
%   ct = matRad_addmovement(ct, ct.motionPeriod, ct.numOfCtScen, amp)
%
% input
%   ct :            ct cube
%   pln:            matRad plan meta information struct
%   dij:            matRad dij struct
%   stf:            matRad steering information struct
%   cst:            matRad cst struct
%   resultGUI:      struct containing optimized fluence vector
% output
%   resultGUI:      structure containing phase dose, RBE weighted dose, etc
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a time sequence for when each bixel is irradiated, the sequence
% follows the backforth spot scanning
timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI);

if isfield(ct, 'motionPeriod')
    motionPeriod = ct.motionPeriod;
else
    error('motionPeriod is not defined in ct structure')
end

% prepare a phase matrix
motion = 'linear'; % the assumed motion type
timeSequence = matRad_makePhaseMatrix(timeSequence, ct.numOfCtScen, motionPeriod, motion);

resultGUI.bioParam = pln.bioParam;

totalPhaseMatrix = resultGUI.w*[0 1 0];%vertcat(timeSequence.phaseMatrix);

for iPhase = 1:ct.numOfCtScen
    
    resultGUI.phaseDose{iPhase} = reshape(dij.physicalDose{iPhase} * totalPhaseMatrix(:,iPhase), dij.dimensions);
    
    if isequal(resultGUI.bioParam.model,'MCN')
        
        resultGUI.phaseAlphaDose{iPhase} = reshape(dij.mAlphaDose{iPhase} * w, dij.dimensions);
        resultGUI.phaseSqrtBetaDose{iPhase} = reshape(dij.mSqrtBetaDose{iPhase} * w, dij.dimensions);
        
    end
end

% dose accumulation
resultGUI = matRad_doseAcc(ct, resultGUI, cst, 'DDM');  %acc Methods: 'EMT' 'DDM'
