function [resultGUI, timeSequence] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, totalPhaseMatrix)
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('totalPhaseMatrix','var')
   % make a time sequence for when each bixel is irradiated, the sequence
   % follows the backforth spot scanning
   timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI);

   % prepare a phase matrix
   motion       = 'linear'; % the assumed motion type
   timeSequence = matRad_makePhaseMatrix(timeSequence, ct.numOfCtScen, ct.motionPeriod, motion);

   resultGUI.bioParam = pln.bioParam;

   % the total phase matrix determines what beamlet will be administered in what ct phase
   totalPhaseMatrix   = vertcat(timeSequence.phaseMatrix);
else
   timeSequence = [];
   
end
% compute all phases
for i = 1:ct.numOfCtScen
    
    tmpResultGUI = matRad_calcCubes(totalPhaseMatrix(:,i),dij,i);
    
    % compute physical dose for physical opt
    if strcmp(pln.bioParam.model,'none')       
        resultGUI.phaseDose{i} = tmpResultGUI.physicalDose;
    % compute RBExD with const RBE
    elseif strcmp(pln.bioParam.model,'constRBE')
        resultGUI.phaseRBExD{i} = tmpResultGUI.RBExD;
    % compute all fields
    elseif any(strcmp(pln.bioParam.model,{'MCN','LEM','WED'}))
        resultGUI.phaseAlphaDose{i}    = tmpResultGUI.alpha .* tmpResultGUI.physicalDose;
        resultGUI.phaseSqrtBetaDose{i} = sqrt(tmpResultGUI.beta) .* tmpResultGUI.physicalDose;
    end
    
end

% accumulation
if strcmp(pln.bioParam.model,'none')       
    
    resultGUI.accPhysicalDose = matRad_doseAcc(ct,resultGUI.phaseDose, cst, 'DDM');

elseif strcmp(pln.bioParam.model,'constRBE')

    resultGUI.accRBExD = matRad_doseAcc(ct,resultGUI.phaseRBExD, cst, 'DDM');

elseif any(strcmp(pln.bioParam.model,{'MCN','LEM','WED'}))

    resultGUI.accAlphaDose    = matRad_doseAcc(ct,resultGUI.phaseAlphaDose, cst, 'DDM');
    resultGUI.accSqrtBetaDose = matRad_doseAcc(ct,resultGUI.phaseSqrtBetaDose, cst, 'DDM');

    % only compute where we have biologically defined tissue
    ix = dij.alphaX~=0; 
    
    resultGUI.accEffect = resultGUI.accAlphaDose + resultGUI.accSqrtBetaDose.^2;
    
    resultGUI.accRBExD     = zeros(ct.cubeDim);
    resultGUI.accRBExD(ix) = ((sqrt(dij.alphaX(ix).^2 + 4 .* dij.betaX(ix) .* resultGUI.accEffect(ix)) - dij.alphaX(ix))./(2.*dij.betaX(ix)));
        
end

