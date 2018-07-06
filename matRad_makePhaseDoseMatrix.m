function phaseMatrix = matRad_makePhaseDoseMatrix(resultGUI, bixelInfo, numOfPhases, motionPeriod)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% call
%   
%
% input
%       
%  
% output
%
% comment:
% 
% References
%
% Silke Ulrich May 2017
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phaseTime = motionPeriod/numOfPhases;

for i = 1:length(bixelInfo)
    bixelInfo(i).phaseMatrix = zeros(length(bixelInfo(i).time),numOfPhases);
end

