function timeSequence = matRad_makePhaseMatrix(timeSequence, numOfPhases, motionPeriod, motion)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   using the time sequence and the ordering of the bixel iradiation, and
%   number of scenarios, makes a phase matrix of size number of bixels *
%   number of scenarios 
%
%
% call
%   timeSequence = matRad_makePhaseMatrix(timeSequence, numOfPhases, motionPeriod, motion)
%
% input
%   timeSequence:   struct containing bixel ordering information and the
%                   time sequence of the spot scanning
%   numOfCtScen:    number of the desired phases
%   motionPeriod:   the extent of a whole breathing cycle (in seconds)
%   motion:         motion scenario: 'linear'(default), 'sampled_period' 
%
% output
%   timeSequence:      phase matrix field added
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

if nargin < 4
    motion = 'linear';
end

phaseTimeDev = .01;
% time of each phase in micro seconds
phaseTime = motionPeriod * 10 ^ 6/numOfPhases;

for i = 1:length(timeSequence)
    
    realTime = phaseTime;
    timeSequence(i).phaseMatrix = zeros(length(timeSequence(i).time),numOfPhases);
    
    iPhase = 1;
    iTime = 1;

    while (iTime <= length(timeSequence(i).time))
        if(timeSequence(i).time(iTime) < realTime)
            
            while(iTime <= length(timeSequence(i).time) && timeSequence(i).time(iTime) < realTime)
                timeSequence(i).phaseMatrix(iTime, iPhase) = 1;
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
    
    % permuatation of phaseMatrix from SS order to STF order
    timeSequence(i).phaseMatrix = timeSequence(i).phaseMatrix(timeSequence(i).orderToSTF,:);
    
    % inserting the fluence in phaseMatrix
    timeSequence(i).phaseMatrix = timeSequence(i).phaseMatrix .* timeSequence(i).w;
    
end

