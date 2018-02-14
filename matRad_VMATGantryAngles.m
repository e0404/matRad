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
%   numApertures
%   numLevels
%   minGantryAngleRes
%   maxApertureAngleSpread



gantryAngleSpacing = pln.minGantryAngleRes; %ideally should be spaced every 2 or 4 degrees; gantry spacing that dij is performed

maxNumApertures = pln.maxApertureAngleSpread/gantryAngleSpacing;

gantryToOptAngleSpacingFactor = floor(maxNumApertures/pln.numApertures);
optGantryAngleSpacing = gantryAngleSpacing*gantryToOptAngleSpacingFactor;
initGantryAngleSpacing = pln.numApertures*optGantryAngleSpacing;

pln.gantryAngles    = 0:gantryAngleSpacing:360;
pln.optGantryAngles = 0:optGantryAngleSpacing:360;
pln.initGantryAngles = (optGantryAngleSpacing*floor(pln.numApertures/2)):initGantryAngleSpacing:361;

if pln.gantryAngles(end) ~= pln.optGantryAngles(end)
    %the final beam is not an optimized beam yet, but it should be for
    %interpolation reasons
    pln.optGantryAngles(numel(pln.optGantryAngles)+1) = pln.gantryAngles(end);
end

if pln.initGantryAngles(1) == 0 && pln.initGantryAngles(end) ~= 360
    pln.initGantryAngles(numel(pln.initGantryAngles)+1) = 360;
end

pln.numOfBeams      = numel(pln.gantryAngles);

pln.couchAngles     = 0*pln.gantryAngles;


pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);



