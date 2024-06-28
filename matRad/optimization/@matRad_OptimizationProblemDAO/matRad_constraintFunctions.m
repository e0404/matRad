function c = matRad_constraintFunctions(optiProb,apertureInfoVec,dij,cst)
% matRad IPOPT callback: constraint function for direct aperture optimization
% 
% call
%   c = matRad_constraintFunctions(optiProb,apertureInfoVec,dij,cst)
%
% input
%   optiProb:       option struct defining the type of optimization
%   apertueInfoVec: aperture info vector
%   dij:            dose influence matrix
%   cst:            matRad cst struct
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

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
if ~isequal(apertureInfoVec,optiProb.apertureInfo.apertureVector)
    optiProb.apertureInfo = optiProb.matRad_daoVec2ApertureInfo(optiProb.apertureInfo,apertureInfoVec);
end
apertureInfo = optiProb.apertureInfo;

% value of constraints for leaves
leftLeafPos  = apertureInfoVec([1:apertureInfo.totalNumOfLeafPairs]+apertureInfo.totalNumOfShapes);
rightLeafPos = apertureInfoVec(1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes:end);
c_dao        = rightLeafPos - leftLeafPos;

% bixel based objective function calculation
c_dos = matRad_constraintFunctions@matRad_OptimizationProblem(optiProb,apertureInfo.bixelWeights,dij,cst);

% concatenate
c = [c_dao; c_dos];
