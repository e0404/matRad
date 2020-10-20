function stf = matRad_generateStfPencilBeam(pln,ct,energyIx)
% function to generate very basic stf file for a single energy
% 
% call
%   stf = matRad_generateStfPencilBeam(pln,EnergyIx)
%
% input
%   pln:                    matRads planning struct
%   energyIx (optional):    index of desired energy of the machine data
%                           warning: only optional if isocenter defined in
%                           pln and beam/couch angles = 0;
%
% output
%   stf:            matRads steering struct
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global matRad_cfg;
matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('matRad: Generating pencil beam ...\n');

load([pln.radiationMode,'_',pln.machine]);

SAD = machine.meta.SAD;
currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:)', machine.meta.LUT_bxWidthminFWHM(2,:)',pln.propStf.bixelWidth);

% get roation matrix for gantry and couch angles unequal to 0
rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles,pln.propStf.couchAngles));

% generate stf parameters
stf.gantryAngle = pln.propStf.gantryAngles;
stf.couchAngle = pln.propStf.couchAngles;
stf.bixelWidth = pln.propStf.bixelWidth;
stf.radiationMode = pln.radiationMode;
stf.SAD = SAD;
stf.isoCenter = pln.propStf.isoCenter;

stf.numOfRays = 1;
stf.numOfBixelsPerRay = 1;
stf.totalNumOfBixels = 1;

stf.sourcePoint_bev = [0,-SAD,0];
stf.sourcePoint = stf.sourcePoint_bev*rotMat_vectors_T;

stf.ray.rayPos_bev = [0,0,0];
stf.ray.targetPoint_bev = [0,SAD,0];
stf.ray.rayPos = stf.ray.rayPos_bev*rotMat_vectors_T;
stf.ray.targetPoint = stf.ray.targetPoint_bev*rotMat_vectors_T;

% automatic energy selection if no energy selected
if nargin < 3
    [~,l{1},rho{1},~,~] = matRad_siddonRayTracer(pln.propStf.isoCenter + pln.multScen.isoShift(1,:), ...
        ct.resolution, stf.sourcePoint, stf.ray.targetPoint, [ct.cube]);
    radDepths = cumsum(l{1} .* rho{1}{1});
    currRadDepths = radDepths(80);
    [~,energyIx] = min(abs([machine.data.peakPos]-currRadDepths));
end

% generate ray
stf.ray.energy = machine.data(energyIx).energy;
stf.ray.focusIx = find(machine.data(energyIx).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
stf.ray.rangeShifter = struct;
stf.ray.rangeShifter.ID = 0;
stf.ray.rangeShifter.eqThickness = 0;
stf.ray.rangeShifter.sourceRashiDistance = 0;

end