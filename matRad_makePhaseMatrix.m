function bixelInfo= matRad_makePhaseMatrix(bixelInfo, numOfPhases, motionPeriod, motion)
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

phaseTimeDev = .01;
phaseTime = motionPeriod * 10 ^ 6/numOfPhases;



for i = 1:length(bixelInfo)
    
    realTime = phaseTime;
    bixelInfo(i).phaseMatrix = zeros(length(bixelInfo(i).time),numOfPhases);
    
    iPhase = 1;
    iTime = 1;

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
            
            switch motion
                case 'linear'
                    ...
                case 'sampled_period'
                    phaseTime = phaseTime + phaseTimeDev * randn(1);
                case 'periodic'
                    % To Do
                otherwise
                    fprintf('The specified motion scenario is not included.')
            end
            
            realTime = realTime + phaseTime;
        end
    end
end

%         phaseTime = motionPeriod * 10 ^ 6/numOfPhases;
%         amp = phaseTime * motionPeriod / 4;
%         periodicTerm =  - phaseTime / 2 + amp * abs(sin(pi * realTime/ motionPeriod));
%             periodicTerm = abs(sin((pi * realTime  * Ix/ motionPeriod) + pi/2))
