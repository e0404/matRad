function [resultGUI, bixelInfo] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad 4D dose calculation
%
% call
%
%
% input
%
%
% output
%
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 201( the matRad development team. 
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

motionPeriod = 5; % a whole breathing motion period (in seconds)
motion = 'linear'; % the assumed motion of chest
numOfPhases = size(ct.cube, 2);

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

%dose accumulation
resultGUI = matRad_doseAcc(ct, resultGUI, cst, 'DDM');  %acc Methods: 'EMT' 'DDM'

%visualisation
matRad_plotPhaseDose_2(ct, cst, pln, resultGUI); %optional kann slice angegeben werden  

% Plot the result in comparison to the static dose
slice = round(pln.isoCenter(1,3)./ct.resolution.z); 
figure 
subplot(2,2,1)
imagesc(resultGUI.RBExD(:,:,slice)),colorbar, colormap(jet); 
title('static dose distribution')
subplot(2,2,2)
imagesc(resultGUI.accRBExD(:,:,slice)),colorbar, colormap(jet); 
title('4D dose distribution')
subplot(2,2,3)
imagesc(resultGUI.RBExD(:,:,slice) - resultGUI.accRBExD(:,:,slice)) ,colorbar, colormap(jet); 
title('Difference')