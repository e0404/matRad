function pln = matRad_VMATGantryAngles(pln,way)
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

if strcmp(way,'new')
    
    pln.gantryAngleSpacing = pln.minGantryAngleRes; %ideally should be spaced every 2 or 4 degrees; gantry spacing that dij is performed
    
    pln.maxNumApertures = pln.maxApertureAngleSpread/pln.gantryAngleSpacing;
    
    pln.gantryToOptAngleSpacingFactor = floor(pln.maxNumApertures/pln.numApertures);
    pln.optGantryAngleSpacing = pln.gantryAngleSpacing*pln.gantryToOptAngleSpacingFactor;
    pln.initGantryAngleSpacing = pln.numApertures*pln.optGantryAngleSpacing;
    
    pln.gantryAngles    = 0:pln.gantryAngleSpacing:360;
    pln.optGantryAngles = 0:pln.optGantryAngleSpacing:360;
    pln.initGantryAngles = (pln.optGantryAngleSpacing*floor(pln.numApertures/2)):pln.initGantryAngleSpacing:361;
    
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
    
elseif strcmp(way,'old')
    
    
end





