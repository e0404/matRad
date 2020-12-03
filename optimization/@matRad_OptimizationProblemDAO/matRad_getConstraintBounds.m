function [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
% matRad IPOPT get constraint bounds function for direct aperture optimization
% 
% call
%   [cl,cu] = matRad_daoGetConstBounds(optiProb,cst)
%
% input
%   optiProb:   option struct defining the type of optimization
%   cst:        matRad cst struct
%
% output
%   cl:         lower bounds on constraints
%   cu:         lower bounds on constraints
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

apertureInfo = optiProb.apertureInfo;

% Initialize bounds
cl_dao = zeros(apertureInfo.totalNumOfLeafPairs,1);
cu_dao = inf*ones(apertureInfo.totalNumOfLeafPairs,1);

% get dosimetric bounds from cst (just like for conv opt) by call to
% superclass method
[cl_dos,cu_dos] = matRad_getConstraintBounds@matRad_OptimizationProblem(optiProb,cst);

% concatenate
cl = [cl_dao; cl_dos];
cu = [cu_dao; cu_dos];
