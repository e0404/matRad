function [stf, pln] = matRad_importDicomSteeringPhotons(pln)
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

stf = struct;
if ~isfield(pln.propStf.collimation,'Fields')
    return
end

% get fields possessing a field weight vector greater than 0
Fields = pln.propStf.collimation.Fields([pln.propStf.collimation.Fields(:).Weight] > 0);

[UniqueComb,ia,ib] = unique( vertcat([Fields(:).GantryAngle], [Fields(:).CouchAngle])','rows');

% return corret angles to pln, because some angle derivations might be 
% only in the control point sequences
pln.propStf.gantryAngles = UniqueComb(:,1)';
pln.propStf.couchAngles  = UniqueComb(:,2)';

stf = struct;
% loop over all fields
for i = 1:size(UniqueComb,1)
    % set necessary steering information 
    stf(i).gantryAngle  = UniqueComb(i,1);
    stf(i).couchAngle   = UniqueComb(i,2);
    stf(i).isoCenter    = pln.propStf.isoCenter(i,:);
    
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
    
    % coordinate transformation with rotation matrix.
    % use transpose matrix because we are working with row vectors
    rotMat_vectors_T = transpose(matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle));


    % Rotated Source point (1st gantry, 2nd couch)
    stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;
    
    % only one ray in center position
    stf(i).ray.rayPos_bev = [0 0 0];
    stf(i).ray.rayPos     = stf(i).ray.rayPos_bev*rotMat_vectors_T;

    % target point is for ray in center position at
    stf(i).ray.targetPoint_bev = [0 stf(i).SAD 0];
    stf(i).ray.targetPoint     = stf(i).ray.targetPoint_bev*rotMat_vectors_T;
    
    % set weight for output field
    stf(i).ray.weight = 1; % weighting incorporated into primary fluence --> finalShape
    %stf(i).ray.SSD    = Fields(ia(i)).SSD;
    stf(i).ray.energy = Fields(ia(i)).Energy;  
    
    ix = (ib == i);
    currFieldSeq = Fields(ix);
    
    % add weighted shapes for the current beam
    finalShape = 0;
    for j = 1:sum(ix)
        finalShape = finalShape + currFieldSeq(j).Weight * currFieldSeq(j).Shape;
    end
    stf(i).ray.shape = finalShape;

end

