classdef matRad_PriorityListConstraint < handle
% matRad_PriorityListConstraint implements a helper class for the priority list
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
        constraint; 
        cstIdx;    % idx of the volume in the cst
    end

    methods 
        function obj = matRad_PriorityListConstraint(constraint,cstIdx)
            %constructor
            obj.constraint = constraint;
            obj.cstIdx = cstIdx;
        end
    end
end