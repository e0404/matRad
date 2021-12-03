function [DistanceMatrix,DistanceVector] = matRad_getDistanceMatrix(seedPoints,dosePoints)
% matRad_getDistanceMatrix gets (seedpoint x dosepoint) matrix of relative
% distances
%
% call
%   [DistanceMatrix,DistanceVector] = getDistanceMatrix(seedPoints,...
%       dosePoints)
%   normally called within matRad_getBrachyDose
%
% input
%   seedPoints:     struct with fields x,y,z
%   dosePoints:     struct with fields x,y,z
%
% output
%   distance matrix:    rows: index of dosepoint 
%                       columns: index of deedpoint
%                       entry: distance of seedpoints and dosepoint in mm
%                       |
%                       | DistanceMatrix.x/y/z:   x/y/z component of 
%                           distance(needed for theta calc)
%                       | DistanceMatrix.dist: eucledian distance
%   distance vector:    column vector of DistanceMatrix.dist entries
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



DistanceMatrix.x = dosePoints.x'*ones(1,length(seedPoints.x)) - ones(length(dosePoints.x),1)*seedPoints.x;
DistanceMatrix.y = dosePoints.y'*ones(1,length(seedPoints.y)) - ones(length(dosePoints.y),1)*seedPoints.y;
DistanceMatrix.z = dosePoints.z'*ones(1,length(seedPoints.z)) - ones(length(dosePoints.z),1)*seedPoints.z;
DistanceMatrix.dist = sqrt(DistanceMatrix.x.^2+DistanceMatrix.y.^2+DistanceMatrix.z.^2);

if nargout == 2
DistanceVector = reshape(DistanceMatrix.dist,[],1);
end

end

