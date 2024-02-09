function dominatedPoints = matRad_paretoCheck(fVals,newfVal)
% matRad helper function that checks if a newly generated point dominates
% another point or is dominated by a point in the previous sets. should
% be called after calculating the first initial points 
%
%
% input
%   fVals:              The so-far determined pareto optimal points 
%   newfVal:            New point to check against fVals
%
% output
%   dominatedPoints:    0 if new point dominated by point in previous set
%                       empty if new point neither dominates a point in old set or is dominated by point in old set
%                       otherwise array refering to the indices of points
%                       in the old set dominated by the newly calculated
%                       point
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dominatedPoints = [];
for i = 1:size(fVals,1)
    if all(newfVal > fVals(i,:)) % newly generated point is not pareto optimal
        dominatedPoints = [dominatedPoints,0];
        break 
    elseif all(newfVal < fVals(i,:)) % point i  is dominated by newly generated point
        dominatedPoints = [dominatedPoints,i];
    end
end