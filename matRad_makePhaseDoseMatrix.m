function bixelInfo= matRad_makePhaseDoseMatrix(bixelInfo, numOfPhases, motionPeriod, motion)
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
% Ahmad Neishabouri June 2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    motion = 'linear';
end

iTime = 1;

switch motion
    case 'linear'        
        phaseTime = motionPeriod * 10 ^ 6/numOfPhases;
    case 'periodic'
           % To Do
    otherwise
        fprintf('The specified motion scenario is not included.')
end

realTime = phaseTime;

for i = 1:length(bixelInfo)
    
    bixelInfo(i).phaseMatrix = zeros(length(bixelInfo(i).time),numOfPhases);
    iPhase = 1;

    while (iTime <= length(bixelInfo(i).time))
        if(bixelInfo(i).time_ordered(iTime) < realTime)
            bixelInfo(i).time_ordered(iTime);
            while(iTime <= length(bixelInfo(i).time) && bixelInfo(i).time_ordered(iTime) < realTime)
                bixelInfo(i).phaseMatrix(iTime, iPhase) = 1;
                iTime = iTime + 1;
            end
        else
            iPhase = iPhase + 1;
            if(iPhase > numOfPhases)
                iPhase = 1;
            end
            
            realTime = realTime + phaseTime;
        end
    end
end


