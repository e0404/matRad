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


if ~isfield(pln.propOpt.VMAToptions,'maxGantryAngleSpacing')
    error('Please define pln.propOpt.maxGantryAngleSpacing.');
end

if ~isfield(pln.propOpt.VMAToptions,'maxDAOGantryAngleSpacing')
    error('Please define pln.propOpt.maxDAOGantryAngleSpacing.');
end

if ~isfield(pln.propOpt.VMAToptions,'maxFMOGantryAngleSpacing')
    error('Please define pln.propOpt.maxFMOGantryAngleSpacing.');
end

% ensure that an integer number of DAO gantry angles will fit in 360 degrees,
% spaced at least as close as maxDAOGantryAngleSpacing
numDAOGantryAngles = ceil(360/pln.propOpt.VMAToptions.maxDAOGantryAngleSpacing);
DAOGantryAngleSpacing = 360/numDAOGantryAngles;
% actual number of gantry angles is numDAOGantryAngles+1;

% ensure that DAOGantryAngleSpacing is an integer multiple of gantryAngleSpacing
numGantryAngles = ceil(numDAOGantryAngles*DAOGantryAngleSpacing/pln.propOpt.VMAToptions.maxGantryAngleSpacing);
gantryAngleSpacing = 360/numGantryAngles;

% ensure that FMOGantryAngleSpacing is an odd integer multiple of DAOGantryAngleSpacing
numApertures = floor(pln.propOpt.VMAToptions.maxFMOGantryAngleSpacing/DAOGantryAngleSpacing);
if mod(numApertures,2) == 0
    numApertures = numApertures-1;
end
FMOGantryAngleSpacing = numApertures*DAOGantryAngleSpacing;
firstFMOGantryAngle = DAOGantryAngleSpacing*floor(numApertures/2);
lastFMOGantryAngle = 360-DAOGantryAngleSpacing*floor(numApertures/2);


pln.propStf.gantryAngles    = 0:gantryAngleSpacing:360;
pln.propStf.DAOGantryAngles = 0:DAOGantryAngleSpacing:360;
pln.propStf.FMOGantryAngles = firstFMOGantryAngle:FMOGantryAngleSpacing:lastFMOGantryAngle;


pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.couchAngles     = 0*pln.propStf.gantryAngles;
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);



