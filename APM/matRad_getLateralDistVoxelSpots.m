function [ vLatDistX,vLatDistZ ] = matRad_getLateralDistVoxelSpots(ix,cubeDim,CTres,isoCenter,stf,pln,vCandidates,SAD)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the lateral distance from a given voxel to the cenetral axis of defined spots
% 
% call
%   [ vLatDistX,vLatDistZ ] = matRad_getLateralDistVoxelSpots(ix,cubeDim,CTres,isoCenter,stf,pln,vCandidates,SAD)
%
% input
%   ix:               voxel index
%   cubeDim:          cube dimensions
%   isoCenter:        iso center
%   stf:              matRad's stf struct
%   pln:              matRad's pln struct
%   vCandidates;      spot index vector
%   SAD;              source to axis distance --> machine.meta.SAD
%
% output
%   vLatDistX:        lateral distance in x from voxel to central ray
%   vLatDistZ:        lateral distance in z from voxel to central ray
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

USE_MEX = false;   

if isunix && ~ismac
   USE_MEX = false; % mex files are available for mac and windows but not for linux
end

% check if mex file is available
if exist('matRad_calcGeoDistsFast_mex','file') == 3
    USE_MEX = false;
end

[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(cubeDim,ix);
xCoordsV = xCoordsV_vox(:) * CTres(1) - isoCenter(1);
yCoordsV = yCoordsV_vox(:) * CTres(2) - isoCenter(2);
zCoordsV = zCoordsV_vox(:) * CTres(3) - isoCenter(3);
coordsV  = [xCoordsV yCoordsV zCoordsV];

maxLateralCutoffDoseCalc = 70;

radDepthIx = 1; 
counter    = 1;
numRay     = pln.multScen.bixelRayLUT(vCandidates);
numBeam    = pln.multScen.bixelBeamLUT(vCandidates);
vLatDistX  = zeros(numel(vCandidates),1); 
vLatDistZ  = zeros(numel(vCandidates),1);


for k = 1:length(vCandidates)
   
 % get rotation matrix
 rotMat_system_T = matRad_getRotationMatrix(pln.propStf.gantryAngles(numBeam(k)),pln.propStf.couchAngles(numBeam(k)));
 % rotate coordinates (1st couch around Y axis, 2nd gantry movement)
 rot_coordsV      = (coordsV*rotMat_system_T);
 rot_coordsV(:,1) = rot_coordsV(:,1) - stf(numBeam(k)).sourcePoint_bev(1);
 rot_coordsV(:,2) = rot_coordsV(:,2) - stf(numBeam(k)).sourcePoint_bev(2);
 rot_coordsV(:,3) = rot_coordsV(:,3) - stf(numBeam(k)).sourcePoint_bev(3);

 if USE_MEX

    [vLatDistX(counter),vLatDistZ(counter)] = matRad_calcGeoDistsFast_mex(rot_coordsV, ...
                                              stf(numBeam(k)).sourcePoint_bev, ...
                                              stf(numBeam(k)).ray(numRay(k)).targetPoint_bev, ...
                                              SAD, ...
                                              radDepthIx, ...
                                              maxLateralCutoffDoseCalc);
 else
    [vLatDistX(counter),vLatDistZ(counter)] = matRad_calcGeoDistsFast(rot_coordsV, ...
                                              stf(numBeam(k)).sourcePoint_bev, ...
                                              stf(numBeam(k)).ray(numRay(k)).targetPoint_bev, ...
                                              SAD, ...
                                              radDepthIx, ...
                                              maxLateralCutoffDoseCalc);
  end % eof MEX            

   counter = counter + 1;
   
end
             
                                                    


end % eof function

