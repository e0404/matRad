function [stf, pln] = matRad_importDicomSteeringPhotons(pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to import a matRad stf struct from dicom RTPLAN data
% 
% call
%   [stf, pln] = matRad_importDicomSteeringPhotons(pln)
%
% input
%   pln:            matRad pln struct with meta information (collimation
%                   data included)
%
% output
%   stf             matRad stf struct
%   pln             matRad pln struct
%
% References
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
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

% get fields possessing a field weight vector greater than 0
Fields = pln.Collimation.Fields([pln.Collimation.Fields(:).Weight] > 0);

[UniqueComb,ia,ib] = unique( vertcat([Fields(:).GantryAngle], [Fields(:).CouchAngle])','rows');

% return corret angles to pln, because some angle derivations might be 
% only in the control point sequences
pln.gantryAngles = UniqueComb(:,1)';
pln.couchAngles = UniqueComb(:,2)';

stf = struct;
% loop over all fields
for i = 1:size(UniqueComb,1)
    % set necessary steering information 
    stf(i).gantryAngle  = UniqueComb(i,1);
    stf(i).couchAngle   = UniqueComb(i,2);
    stf(i).isoCenter    = pln.isoCenter;
    
    % bixelWidth = 'field' as keyword for whole field dose calc
    stf(i).bixelWidth    = 'field';
    stf(i).radiationMode = 'photons';
    
    % only one bixel per ray and one ray for photon dose calc based on
    % fields
    stf(i).numOfBixelsPerRay = 1;
    stf(i).numOfRays         = 1;
    stf(i).totalNumOfBixels  = stf(i).numOfRays;
    stf(i).SAD               = Fields(ia(i)).SAD;
    stf(i).sourcePoint_bev   = [0 -stf(i).SAD 0];
    
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient; use transpose matrices because we are
    % working with row vectors
    % Rotation around Z axis (gantry)
    rotMx_XY_T = [ cosd(stf(i).gantryAngle) sind(stf(i).gantryAngle) 0;
                  -sind(stf(i).gantryAngle) cosd(stf(i).gantryAngle) 0;
                                          0                        0 1];
    % Rotation around Y axis (couch)
    rotMx_XZ_T = [cosd(stf(i).couchAngle) 0 -sind(stf(i).couchAngle);
                                        0 1                        0;
                  sind(stf(i).couchAngle) 0  cosd(stf(i).couchAngle)];

    % Rotated Source point (1st gantry, 2nd couch)
    stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMx_XY_T*rotMx_XZ_T;
    
    % only one ray in center position
    stf(i).ray.rayPos_bev = [0 0 0];
    stf(i).ray.rayPos     = stf(i).ray.rayPos_bev*rotMx_XY_T*rotMx_XZ_T;

    % target point is for ray in center position at
    stf(i).ray.targetPoint_bev = [0 stf(i).SAD 0];
    stf(i).ray.targetPoint     = stf(i).ray.targetPoint_bev*rotMx_XY_T*rotMx_XZ_T;
    
    % set weight for output field
    stf(i).ray.weight = Fields(ia(i)).FinalCumWeight;
    stf(i).ray.SSD    = Fields(ia(i)).SSD;
    stf(i).ray.energy = Fields(ia(i)).Energy;  
    
    ix = (ib == i);
    currFieldSeq = Fields(ix);
    
    % add weighted shapes for the current beam
    partialWeight   = [currFieldSeq(:).Weight]./stf(i).ray.weight;
    finalShape = 0;
    for j = 1:sum(ix)
        finalShape = finalShape + partialWeight(j) * currFieldSeq(j).Shape;
    end
    stf(i).ray.shape = finalShape;
end

end

