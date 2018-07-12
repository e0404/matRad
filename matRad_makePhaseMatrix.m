function bixelInfo = matRad_makePhaseMatrix(bixelInfo, numOfPhases, motionPeriod, motion)
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
        if(bixelInfo(i).time(iTime) < realTime)
            
            bixelInfo(i).time(iTime);
            
            while(iTime <= length(bixelInfo(i).time) && bixelInfo(i).time(iTime) < realTime)
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
    
    bixelInfo(i).phaseMatSTF = bixelInfo(i).phaseMatrix(bixelInfo(i).orderToSTF,:);
    bixelInfo(i).phaseMatSTF = bixelInfo(i).phaseMatSTF .* bixelInfo(i).w;
    
end

bixelInfo(1).phaseMatSTF_total = vertcat(bixelInfo.phaseMatSTF);