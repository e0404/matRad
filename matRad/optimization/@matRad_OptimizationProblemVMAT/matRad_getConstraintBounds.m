function [cl, cu] = matRad_getConstraintBounds(optiProb, cst)
% matRad IPOPT get constraint bounds function for VMAT
%
% call
%   [cl,cu] = matRad_daoGetConstBounds(cst,apertureInfo,type)
%
% input
%   cst:            matRad cst struct
%   apertureInfo:   aperture info struct
%   options:        option struct defining the type of optimization
%
% output
%   cl: lower bounds on constraints
%   cu: lower bounds on constraints
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

% get dosimetric bounds from cst by call to DAO superclass method
[cl_dos_dao, cu_dos_dao] = matRad_getConstraintBounds@matRad_OptimizationProblemDAO(optiProb, cst);

optInd = find([apertureInfo.arc.beam.isDAOBeam]);

% numOfActiveLeafPairs should be independent of the beam, due to using the
% union of all ray positions in the stf
nLeafPairs = apertureInfo.beam(1).numOfActiveLeafPairs;
leafSpeedLim = apertureInfo.constraints.leafSpeed;
muRateLim = apertureInfo.constraints.monitorUnitRate;

% Convert from cm/deg when checking constraints; cannot do it at this stage
% since gantry rotation speed is not hard-coded
if apertureInfo.continuousAperture
    nLeafSpeed = 2 * apertureInfo.arc.numLeafSpeedConstraint * nLeafPairs;
else
    nLeafSpeed = 2 * (numel(optInd) - 1) * nLeafPairs;
end
cl_lfspd = leafSpeedLim(1) * ones(nLeafSpeed, 1); % Minimum leaf travel speed (mm/s)
cu_lfspd = leafSpeedLim(2) * ones(nLeafSpeed, 1); % Maximum leaf travel speed (mm/s)

cl_dosrt = muRateLim(1) * ones(numel(optInd), 1); % Minimum MU/sec
cu_dosrt = muRateLim(2) * ones(numel(optInd), 1); % Maximum MU/sec

% concatenate
cl = [cl_dos_dao; cl_lfspd; cl_dosrt];
cu = [cu_dos_dao; cu_lfspd; cu_dosrt];
