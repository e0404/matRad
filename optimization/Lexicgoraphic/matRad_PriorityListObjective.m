classdef matRad_PriorityListObjective < handle
% matRad_PriorityListObjective implements a helper class for a priority list
% object
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        objective; 
        goalValue; % aspired value for the objective
        cstIdx;    % idx of the volume in the cst
        achievedValue = -1; %not assigned yet
        achievedValue2 = -1; %not assigned yet
    end

    properties (Access = private)
        VOIIdx;    % idx of the objective for the VOI   
    end

    methods 
        function obj = matRad_PriorityListObjective(objective,goal,cstIdx)
            %constructor
            obj.objective = objective;
            obj.goalValue = goal;
            obj.cstIdx = cstIdx;
        end

        function setVOIIdx(obj,idx)
            obj.VOIIdx = idx;
        end
        
        function idx = getVOIIdx(obj)
            idx = obj.VOIIdx;
        end
    end
end