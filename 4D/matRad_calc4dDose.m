function [resultGUI, bixelInfo] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI)
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
bixelInfo = matRad_makeBixelTimeSeq(stf, resultGUI);


motion = 'linear'; % the assumed motion of chest
numOfPhases = size(ct.cube, 2);

if isfield(ct, 'motionPeriod')
    motionPeriod = ct.motionPeriod;
else
    disp('motionPeriod is not defined in ct structure')
    % TODO: make the above as matRad error message
end

% prepare a phase matrix in which place each bixel dose in it's phase
bixelInfo = matRad_makePhaseMatrix(bixelInfo, numOfPhases, motionPeriod, motion);

resultGUI.bioParam = pln.bioParam;

for iPhase = 1:numOfPhases
    
    w = bixelInfo(1).totalPhaseMatrix(:,iPhase);
    
    resultGUI.phaseDose{iPhase} = reshape(dij.physicalDose{iPhase} * w, dij.dimensions);
    
    if isequal(resultGUI.bioParam.model,'MCN')
        
        resultGUI.phaseAlphaDose{p} = reshape(dij.mAlphaDose{p} * w, dij.dimensions);
        resultGUI.phaseSqrtBetaDose{p} = reshape(dij.mSqrtBetaDose{p} * w, dij.dimensions);
        
    end
end

% dose accumulation
resultGUI = matRad_doseAcc(ct, resultGUI, cst, 'DDM');  %acc Methods: 'EMT' 'DDM'
