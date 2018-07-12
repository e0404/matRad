function [resultGUI, bixelInfo] = matRad_calc4dDoseNish(ct, pln, dij, stf, cst, resultGUI)
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
    
    w = bixelInfo(1).phaseMatSTF_total(:,iPhase);
    
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