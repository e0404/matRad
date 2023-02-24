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

%should have already defined fields:
%   maxGantryAngleSpacing
%   maxDAOGantryAngleSpacing
%   maxFMOGantryAngleSpacing
% if any of these aren't defined, use corresponding DEFAULT value

DEFAULT_GANTRY_SPACING     = 5;
DEFAULT_DAO_GANTRY_SPACING = 10;
DEFAULT_FMO_GANTRY_SPACING = 45;

if ~isfield(pln.propOpt.VMAToptions,'maxGantryAngleSpacing')
    matRad_dispToConsole(['matRad_VMATGantryAngles: Using default value pln.propOpt.maxGantryAngleSpacing = ' num2str(DEFAULT_GANTRY_SPACING) ' degree. \n'],[],'warning');
    pln.propOpt.VMAToptions.maxGantryAngleSpacing    = DEFAULT_GANTRY_SPACING;           % upper bound for gantry angle spacing for dose calculation
end

if ~isfield(pln.propOpt.VMAToptions,'maxDAOGantryAngleSpacing')
    matRad_dispToConsole(['matRad_VMATGantryAngles: Using default value pln.propOpt.maxDAOGantryAngleSpacing= ' num2str(DEFAULT_DAO_GANTRY_SPACING) ' degree. \n'],[],'warning');
    pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing = DEFAULT_DAO_GANTRY_SPACING;     % upper bound for gantry angle spacing for DAO
end

if ~isfield(pln.propOpt.VMAToptions,'maxFMOGantryAngleSpacing')
    matRad_dispToConsole(['matRad_VMATGantryAngles: Using default value pln.propOpt.maxFMOGantryAngleSpacing= ' num2str(DEFAULT_FMO_GANTRY_SPACING) ' degree. \n'],[],'warning');
    pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing = DEFAULT_FMO_GANTRY_SPACING;     % upper bound for gantry angle spacing for FMO
end

angularRange = abs(pln.propOpt.VMAToptions.finishingAngle-pln.propOpt.VMAToptions.startingAngle);

if pln.propOpt.VMAToptions.continuousAperture
    
    % angularRange = gantryAngleSpacing*numGantryAngles
    % ensure that gantryAngleSpacing < maxGantryAngleSpacing (as close as
    % possible)
    numGantryAngles = ceil(angularRange./pln.propOpt.VMAToptions.maxGantryAngleSpacing);
    gantryAngleSpacing = angularRange./numGantryAngles;
    
    % (numDAOGantryAngles-1)*DAOGantryAngleSpacing = (numGantryAngles-1)*gantryAngleSpacing
    % where
    % ensure that DAOGantryAngleSpacing < maxDAOGantryAngleSpacing (as close as
    % possible)
    numDAOGantryAngles = ceil((numGantryAngles-1).*gantryAngleSpacing./pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing)+1;
    % now ensure that numGantryAngles-1 is a multiple of numDAOGantryAngles-1 so
    % that they align
    numGantryAngles = (numDAOGantryAngles-1).*ceil((numGantryAngles-1)./(numDAOGantryAngles-1))+1;
    gantryAngleSpacing = angularRange./numGantryAngles;
    DAOGantryAngleSpacing = (angularRange-gantryAngleSpacing)/(numDAOGantryAngles-1);
    
    % first and last gantry angles are in centre of arc
    firstGantryAngle = pln.propOpt.VMAToptions.startingAngle+gantryAngleSpacing/2;
    lastGantryAngle = pln.propOpt.VMAToptions.finishingAngle-gantryAngleSpacing/2;
    
else
    
    % ensure that an integer number of DAO gantry angles will fit in angularRange degrees,
    % spaced at least as close as maxDAOGantryAngleSpacing
    numDAOGantryAngles = ceil(angularRange/pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing);
    DAOGantryAngleSpacing = angularRange/numDAOGantryAngles;
    % actual number of gantry angles is numDAOGantryAngles+1;
    % ensure that DAOGantryAngleSpacing is an integer multiple of gantryAngleSpacing
    numGantryAngles = ceil(numDAOGantryAngles*DAOGantryAngleSpacing/pln.propOpt.VMAToptions.maxGantryAngleSpacing);
    gantryAngleSpacing = angularRange/numGantryAngles;
    
    % first and last gantry angles are at beginning and end of arc
    firstGantryAngle = pln.propOpt.VMAToptions.startingAngle;
    lastGantryAngle = pln.propOpt.VMAToptions.finishingAngle;
end

% ensure that FMOGantryAngleSpacing is an odd integer multiple of DAOGantryAngleSpacing
numApertures = floor(pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing/DAOGantryAngleSpacing);
if mod(numApertures,2) == 0
    numApertures = numApertures-1;
end
FMOGantryAngleSpacing = numApertures*DAOGantryAngleSpacing;

firstFMOGantryAngle = firstGantryAngle+DAOGantryAngleSpacing*floor(numApertures/2);
lastFMOGantryAngle = lastGantryAngle-DAOGantryAngleSpacing*floor(numApertures/2);

% define angles
pln.propStf.gantryAngles    = firstGantryAngle:gantryAngleSpacing:lastGantryAngle;
pln.propStf.DAOGantryAngles = firstGantryAngle:DAOGantryAngleSpacing:lastGantryAngle;
pln.propStf.FMOGantryAngles = firstFMOGantryAngle:FMOGantryAngleSpacing:lastFMOGantryAngle;

% everything else
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.couchAngles     = 0*pln.propStf.gantryAngles;
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);


%[ixMember,~] = ismember(pln.propStf.gantryAngles,pln.propStf.DAOGantryAngles);
%if sum(ixMember) <= 2

% This should give the expected bevaviour according to the warning message.
if numel(pln.propStf.FMOGantryAngles) <= 2
   matRad_dispToConsole('matRad_VMATGantryAngles: Two or less beam angles will be used for FMO. \n',[],'warning') 
end
