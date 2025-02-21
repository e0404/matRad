function obj = matRad_importDicomSteeringPhotons(obj)
% matRad function to import a matRad stf struct from dicom RTPLAN data
% 
% In your object, there must be a property that contains matRad pln 
% structure with meta information (collimation data included)
% 
% Output - matRad stf and pln structures.
%
% call
%   obj = matRad_importDicomSteeringPhotons(obj)
%
%
% References
%   -
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj.stf = struct;
if ~isfield(obj.pln.propStf.collimation,'Fields')
    return
end

% get fields possessing a field weight vector greater than 0
Fields = obj.pln.propStf.collimation.Fields([obj.pln.propStf.collimation.Fields(:).Weight] > 0);

[UniqueComb,ia,ib] = unique( vertcat([Fields(:).GantryAngle], [Fields(:).CouchAngle])','rows');

% return corret angles to pln, because some angle derivations might be 
% only in the control point sequences
obj.pln.propStf.gantryAngles = UniqueComb(:,1)';
obj.pln.propStf.couchAngles  = UniqueComb(:,2)';

obj.stf = struct;
% loop over all fields
for i = 1:size(UniqueComb,1)
    % set necessary steering information 
    obj.stf(i).gantryAngle  = UniqueComb(i,1);
    obj.stf(i).couchAngle   = UniqueComb(i,2);
    obj.stf(i).isoCenter    = obj.pln.propStf.isoCenter(i,:);
    
    % bixelWidth = 'field' as keyword for whole field dose calc
    obj.stf(i).bixelWidth    = 'field';
    obj.stf(i).radiationMode = 'photons';
    
    % only one bixel per ray and one ray for photon dose calc based on
    % fields
    obj.stf(i).numOfBixelsPerRay = 1;
    obj.stf(i).numOfRays         = 1;
    obj.stf(i).totalNumOfBixels  = obj.stf(i).numOfRays;
    obj.stf(i).SAD               = Fields(ia(i)).SAD;
    obj.stf(i).sourcePoint_bev   = [0 -obj.stf(i).SAD 0];
    
    % coordinate transformation with rotation matrix.
    % use transpose matrix because we are working with row vectors
    rotMat_vectors_T = transpose(matRad_getRotationMatrix(obj.stf(i).gantryAngle,obj.stf(i).couchAngle));


    % Rotated Source point (1st gantry, 2nd couch)
    obj.stf(i).sourcePoint = obj.stf(i).sourcePoint_bev*rotMat_vectors_T;
    
    % only one ray in center position
    obj.stf(i).ray.rayPos_bev = [0 0 0];
    obj.stf(i).ray.rayPos     = obj.stf(i).ray.rayPos_bev*rotMat_vectors_T;

    % target point is for ray in center position at
    obj.stf(i).ray.targetPoint_bev = [0 obj.stf(i).SAD 0];
    obj.stf(i).ray.targetPoint     = obj.stf(i).ray.targetPoint_bev*rotMat_vectors_T;
    
    % set weight for output field
    obj.stf(i).ray.weight = 1; % weighting incorporated into primary fluence --> finalShape
    %obj.stf(i).ray.SSD    = Fields(ia(i)).SSD;
    obj.stf(i).ray.energy = Fields(ia(i)).Energy;  
    
    ix = (ib == i);
    currFieldSeq = Fields(ix);
    
    % add weighted shapes for the current beam
    finalShape = 0;
    for j = 1:sum(ix)
        finalShape = finalShape + currFieldSeq(j).Weight * currFieldSeq(j).Shape;
    end
    obj.stf(i).ray.shape = finalShape;

end

