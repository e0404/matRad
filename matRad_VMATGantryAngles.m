function pln = matRad_VMATGantryAngles(pln,cst,ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad determine gantry angles for VMAT
% 
% call
%   matRad_VMATGantryAngles(pln)
%
% input
%   pln:                matRad plan meta information struct
%
% output
%   pln:                matRad plan meta information struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

USE_DEAULT_VALUES          = false;

DEFAULT_GANTRY_SPACING     = 5;
DEFAULT_DAO_GANTRY_SPACING = 10;
DEFAULT_FMO_GANTRY_SPACING = 45;

if ~isfield(pln.propOpt.VMAToptions,'maxGantryAngleSpacing')
    USE_DEAULT_VALUES = true;
    matRad_dispToConsole(['matRad_VMATGantryAngles: Using default value pln.propOpt.maxGantryAngleSpacing = ' num2str(DEFAULT_GANTRY_SPACING) ' degree. \n'],[],'warning');
end

if ~isfield(pln.propOpt.VMAToptions,'maxDAOGantryAngleSpacing')
    USE_DEAULT_VALUES = true;
     matRad_dispToConsole(['matRad_VMATGantryAngles: Using default value pln.propOpt.maxDAOGantryAngleSpacing= ' num2str(DEFAULT_DAO_GANTRY_SPACING) ' degree. \n'],[],'warning');
end

if ~isfield(pln.propOpt.VMAToptions,'maxFMOGantryAngleSpacing')
    USE_DEAULT_VALUES = true;
     matRad_dispToConsole(['matRad_VMATGantryAngles: Using default value pln.propOpt.maxFMOGantryAngleSpacing= ' num2str(DEFAULT_FMO_GANTRY_SPACING) ' degree. \n'],[],'warning');
end

if USE_DEAULT_VALUES
    pln.propOpt.VMAToptions.maxGantryAngleSpacing    = DEFAULT_GANTRY_SPACING;           % upper bound for gantry angle spacing for dose calculation
    pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = DEFAULT_DAO_GANTRY_SPACING;     % upper bound for gantry angle spacing for DAO
    pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = DEFAULT_FMO_GANTRY_SPACING;     % upper bound for gantry angle spacing for FMO
end

% ensure that an integer number of DAO gantry angles will fit in 360 degrees,
% spaced at least as close as maxDAOGantryAngleSpacing
numDAOGantryAngles    = ceil(360/pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing);
DAOGantryAngleSpacing = 360/numDAOGantryAngles;
% actual number of gantry angles is numDAOGantryAngles+1;

% ensure that DAOGantryAngleSpacing is an integer multiple of gantryAngleSpacing
numGantryAngles    = ceil(numDAOGantryAngles*DAOGantryAngleSpacing/pln.propOpt.VMAToptions.maxGantryAngleSpacing);
gantryAngleSpacing = 360/numGantryAngles;

% ensure that FMOGantryAngleSpacing is an odd integer multiple of DAOGantryAngleSpacing
numApertures = floor(pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing/DAOGantryAngleSpacing);

if mod(numApertures,2) == 0
    numApertures = numApertures-1;
end

FMOGantryAngleSpacing       = numApertures*DAOGantryAngleSpacing;
firstFMOGantryAngle         = DAOGantryAngleSpacing*floor(numApertures/2);
lastFMOGantryAngle          = 360-DAOGantryAngleSpacing*floor(numApertures/2);

pln.propStf.gantryAngles    = 0:gantryAngleSpacing:360;
pln.propStf.DAOGantryAngles = 0:DAOGantryAngleSpacing:360;
pln.propStf.FMOGantryAngles = firstFMOGantryAngle:FMOGantryAngleSpacing:lastFMOGantryAngle;

pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.couchAngles     = 0*pln.propStf.gantryAngles;
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

[ixMember,~] = ismember(pln.propStf.gantryAngles,pln.propStf.DAOGantryAngles);
if sum(ixMember) <= 2
   matRad_dispToConsole('matRad_VMATGantryAngles: Two or less beam angles will be used for FMO. \n',[],'warning') 
end






