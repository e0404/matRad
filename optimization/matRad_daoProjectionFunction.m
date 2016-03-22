function [apertureInfoVec, isConstrActive] = matRad_daoProjectionFunction(apertureInfoVec,limMx,totalNumOfShapes)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to project an iterate during direct aperture optimization
% back onto the feasible set
%
% call
%   [apertureInfoVec, isConstrActive] = matRad_daoProjectionFunction(apertureInfoVec,limMx,totalNumOfShapes)
%
% input
%   apertureInfoVec:  aperture weights and leaf positions stored as vector
%   limMx:            matrix containing upper and lower bounds for apertureInfoVec
%   totalNumOfShapes: total number of aperture shapes
%
% output
%   apertureInfoVec:  vector projected onto feasible set
%   isConstrActive:   boolian vector indicating whether a certain
%                     aperture weight / leaf position is at its boundary
%
% References
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

% project to feasible set
isConstrActive = apertureInfoVec<=limMx(:,1);
apertureInfoVec(apertureInfoVec<=limMx(:,1)) = limMx(apertureInfoVec<=limMx(:,1),1);

isConstrActive = apertureInfoVec>=limMx(:,2) | isConstrActive;
apertureInfoVec(apertureInfoVec>=limMx(:,2)) = limMx(apertureInfoVec>=limMx(:,2),2);

% correct overlapping leaves
totalNumOfLeafPairs = (numel(apertureInfoVec)-totalNumOfShapes)/2;
leftLeafPos = apertureInfoVec(totalNumOfShapes+[1:totalNumOfLeafPairs]);
rightLeafPos = apertureInfoVec(totalNumOfShapes+totalNumOfLeafPairs+[1:totalNumOfLeafPairs]);

overlapping = leftLeafPos >= rightLeafPos;

leftLeafPos(overlapping)  = mean([leftLeafPos(overlapping) rightLeafPos(overlapping)],2);
rightLeafPos(overlapping) = leftLeafPos(overlapping);

apertureInfoVec(totalNumOfShapes+[1:2*totalNumOfLeafPairs]) = [leftLeafPos;rightLeafPos];

overlapping = [zeros(totalNumOfShapes,1);overlapping;overlapping];
isConstrActive = isConstrActive | overlapping;

end

