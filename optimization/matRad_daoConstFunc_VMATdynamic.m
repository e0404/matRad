function c = matRad_daoConstFunc_VMATdynamic(apertureInfoVec,dij,cst,options,daoVec2ApertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: constraint function for direct aperture optimization
%
% call
%   c = matRad_daoObjFunc(apertueInfoVec,dij,cst)
%
% input
%   apertueInfoVec: aperture info vector
%   apertureInfo:   aperture info struct
%   dij:            dose influence matrix
%   cst:            matRad cst struct
%   options: option struct defining the type of optimization
%
% output
%   c:              value of constraints
%
% Reference
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
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

% read in the global apertureInfo and apertureVector variables
global matRad_global_apertureInfo;
% update apertureInfo from the global variable
apertureInfo = matRad_global_apertureInfo;

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
    matRad_global_apertureInfo = apertureInfo;
end

% value of constraints for leaves
leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2);
c_dao        = rightLeafPos - leftLeafPos;

% bixel based objective function calculation
c_dos = matRad_constFuncWrapper(apertureInfo.bixelWeights,dij,cst,options);


% values of times spent in an arc surrounding the optimized angles (full
% arc/dose influence arc)
timeOptBorderAngles = apertureInfoVec((1+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2):end);
timeDoseBorderAngles = timeOptBorderAngles.*[apertureInfo.beam([apertureInfo.beam.optimizeBeam]).timeFacCurr]';

leftLeafDiff = diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);
rightLeafDiff = diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,[]),1,2);

leftLeafDiff = reshape(leftLeafDiff(repmat([apertureInfo.beam.optimizeBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);
rightLeafDiff = reshape(rightLeafDiff(repmat([apertureInfo.beam.optimizeBeam],apertureInfo.beam(1).numOfActiveLeafPairs,1)),apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes);

% values of leaf speeds in an arc surrounding the optimized angles
c_lfspd = reshape(abs([leftLeafDiff rightLeafDiff])./ ...
    repmat(timeDoseBorderAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeDoseBorderAngles),1);

% values of doserate (MU/sec) in an arc surrounding the optimized angles
weights = apertureInfoVec(1:apertureInfo.totalNumOfShapes)./apertureInfo.jacobiScale;
c_dosrt = apertureInfo.weightToMU.*weights./timeDoseBorderAngles;

% concatenate
c = [c_dao; c_lfspd; c_dosrt; c_dos];


