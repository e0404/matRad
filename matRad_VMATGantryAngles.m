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
%   minGantryAngleRes
%   maxApertureAngleSpread

if ~isfield(pln.propOpt.VMAToptions,'minGantryAngleRes')
    error('Please define pln.propOpt.minGantryAngleRes.');
end

if ~isfield(pln.propOpt.VMAToptions,'maxApertureAngleSpread')
    error('Please define pln.propOpt.maxApertureAngleSpread.');
end

if ~isfield(pln.propOpt,'numApertures')
    error('Please define pln.propOpt.numApertures.');
end

if mod(pln.propOpt.VMAToptions.maxApertureAngleSpread,pln.propOpt.VMAToptions.minGantryAngleRes) ~= 0
    error('pln.propOpt.maxApertureAngleSpread should be a multiple of pln.propOpt.minGantryAngleRes.');
end

if mod(pln.propOpt.numApertures,2) ~= 1
    error('pln.propOpt.numApertures should be odd.');
end


gantryAngleSpacing = pln.propOpt.VMAToptions.minGantryAngleRes;

maxNumApertures = pln.propOpt.VMAToptions.maxApertureAngleSpread/gantryAngleSpacing;

if pln.propOpt.numApertures > maxNumApertures
    error('With current settings, pln.propOpt.numApertures should be less than %d.',maxNumApertures);
end

gantryToOptAngleSpacingFactor = floor(maxNumApertures/pln.propOpt.numApertures);
optGantryAngleSpacing = gantryAngleSpacing*gantryToOptAngleSpacingFactor;
initGantryAngleSpacing = pln.propOpt.numApertures*optGantryAngleSpacing;

pln.propStf.gantryAngles    = 0:gantryAngleSpacing:360;
pln.propStf.optGantryAngles = 0:optGantryAngleSpacing:360;
pln.propStf.initGantryAngles = (optGantryAngleSpacing*floor(pln.propOpt.numApertures/2)):initGantryAngleSpacing:361;

if pln.propStf.gantryAngles(end) ~= pln.propStf.optGantryAngles(end)
    %the final beam is not an optimized beam yet, but it should be for
    %interpolation reasons
    pln.propStf.optGantryAngles(numel(pln.optGantryAngles)+1) = pln.propStf.gantryAngles(end);
end

if pln.propStf.initGantryAngles(1) == 0 && pln.propStf.initGantryAngles(end) ~= 360
    pln.propStf.initGantryAngles(numel(pln.propStf.initGantryAngles)+1) = 360;
end


pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);

pln.propStf.couchAngles     = 0*pln.propStf.gantryAngles;

pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);



