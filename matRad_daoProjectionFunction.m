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
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
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

