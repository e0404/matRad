function [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
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
[cl_dos_dao,cu_dos_dao] = matRad_getConstraintBounds@matRad_OptimizationProblemDAO(optiProb,cst);

optInd = find([apertureInfo.propVMAT.beam.DAOBeam]);


if apertureInfo.continuousAperture
    cl_lfspd = apertureInfo.propVMAT.constraints.leafSpeed(1)*ones(2*apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1); %Minimum leaf travel speed (mm/s)
    cu_lfspd = apertureInfo.propVMAT.constraints.leafSpeed(2)*ones(2*apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1); %Maximum leaf travel speed (mm/s)
    %apertureInfo.beam(i).numOfActiveLeafPairs should be independent of i, due to using the union of all ray positions in the stf
    %Convert from cm/deg when checking constraints; cannot do it at this stage since gantry rotation speed is not hard-coded
else
    
    cl_lfspd = apertureInfo.propVMAT.constraints.leafSpeed(1)*ones(2*(numel(optInd)-1)*apertureInfo.beam(1).numOfActiveLeafPairs,1); %Minimum leaf travel speed (mm/s)
    cu_lfspd = apertureInfo.propVMAT.constraints.leafSpeed(2)*ones(2*(numel(optInd)-1)*apertureInfo.beam(1).numOfActiveLeafPairs,1); %Maximum leaf travel speed (mm/s)
    %apertureInfo.beam(i).numOfActiveLeafPairs should be independent of i, due to using the union of all ray positions in the stf
    %Convert from cm/deg when checking constraints; cannot do it at this stage since gantry rotation speed is not hard-coded
end
cl_dosrt = apertureInfo.propVMAT.constraints.monitorUnitRate(1)*ones(numel(optInd),1); %Minimum MU/sec
cu_dosrt = apertureInfo.propVMAT.constraints.monitorUnitRate(2)*ones(numel(optInd),1); %Maximum MU/sec

% concatenate
cl = [cl_dos_dao; cl_lfspd; cl_dosrt];
cu = [cu_dos_dao; cu_lfspd; cu_dosrt];


