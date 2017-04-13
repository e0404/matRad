function apertureInfo = matRad_maxLeafSpeed(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of maximum leaf speed
% 
% call
%   apertureInfo = matRad_maxLeafSpeed(apertureInfo)
%
% input
%   apertureInfo:   aperture info struct
%
% output
%   apertureInfo:   aperture info struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


apertureInfoVec = apertureInfo.apertureVector;

% value of constraints for leaves
leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);

% values of time differences of optimized gantry angles
timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);

i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
j(1) = [];
j(end) = [];

timeFac = [apertureInfo.beam.timeFac]';
timeFac(timeFac == 0) = [];

timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
timeBNOptAngles = timeFacMatrix*timeOptBorderAngles;

% values of leaf speeds of optimized gantry angles
leftLeafSpeed = abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))./repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,1);
rightLeafSpeed = abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2))./repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,1);

% values of max leaf speeds
leftMaxLeafSpeed = max(leftLeafSpeed,[],1);
rightMaxLeafSpeed = max(rightLeafSpeed,[],1);
maxLeafSpeed = max([leftMaxLeafSpeed; rightMaxLeafSpeed],[],1);

% enter into apertureInfo
l = 1;
maxMaxLeafSpeed = 0;
for i = 1:size(apertureInfo.beam,2)
    if apertureInfo.beam(i).optimizeBeam
        
        if l == 1
            apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(l);
        elseif l == apertureInfo.totalNumOfShapes
            apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(l-1);
        else
            apertureInfo.beam(i).maxLeafSpeed = maxLeafSpeed(l-1)*apertureInfo.beam(i).timeFac(1)+maxLeafSpeed(l)*apertureInfo.beam(i).timeFac(2);
        end
        
        
        if l < apertureInfo.totalNumOfShapes && maxLeafSpeed(l) >= maxMaxLeafSpeed
            maxMaxLeafSpeed = maxLeafSpeed(l);
        end
        l = l+1;
    end
end

apertureInfo.maxLeafSpeed = maxMaxLeafSpeed;

