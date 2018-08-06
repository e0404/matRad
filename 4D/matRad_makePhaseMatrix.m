function bixelInfo = matRad_makePhaseMatrix(bixelInfo, numOfPhases, motionPeriod, motion)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   using the time sequence and the ordering of the bixel iradiation, and
%   number of scenarios, makes a phase matrix of size number of bixels *
%   number of scenarios 
%
%
% call
%   bixelInfo = matRad_makePhaseMatrix(bixelInfo, numOfPhases, motionPeriod, motion)
%
% input
%   bixelInfo:      struct containing bixel ordering information and the
%                   time sequence of the spot scanning
%   numOfCtScen:    number of the desired phases
%   motionPeriod:   the extent of a whole breathing cycle (in seconds)
%   motion:         motion scenario: 'linear'(default), 'sampled_period' 
%
% output
%   bixelInfo:      phase matrix field added
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



for i = 1:length(bixelInfo)
    
    
    realTime = phaseTime;
    bixelInfo(i).phaseMatrix = zeros(length(bixelInfo(i).time),numOfPhases);
    
    iPhase = 1;
    iTime = 1;

    while (iTime <= length(bixelInfo(i).time))
        if(bixelInfo(i).time(iTime) < realTime)
            
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
    
    % permuatation of phaseMatrix from SS order to STF order
    bixelInfo(i).phaseMatrix = bixelInfo(i).phaseMatrix(bixelInfo(i).orderToSTF,:);
    % inserting the fluence in phaseMatrix
    bixelInfo(i).phaseMatrix = bixelInfo(i).phaseMatrix .* bixelInfo(i).w;
    
end
% storing all phase matrices in one total matrix
bixelInfo(1).totalPhaseMatrix = vertcat(bixelInfo.phaseMatrix);
