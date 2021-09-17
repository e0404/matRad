function [radDepthV, radDepthCube] = matRad_rayTracing(stfElement,ct,V,rot_coordsV,lateralCutoff)
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
%
% call
%   [radDepthV, radDepthCube] = matRad_rayTracing(stf,ct,V,rot_coordsV,lateralCutoff)
%
% input
%   stfElement:    matRad steering information struct of single(!) beam
%   ct:            ct cube
%   V:             linear voxel indices e.g. of voxels inside patient.
%   rot_coordsV    coordinates in beams eye view inside the patient
%   lateralCutoff: lateral cut off used for ray tracing (optional)

%
% output
%   radDepthV:      radiological depth inside the patient
%   radDepthCube:   radiological depth in whole ct
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S1120179711001359
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
matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispDeprecationWarning('Calls to matRad_rayTracing will be replaced by usage of the matRad_RayTracer class in the future!');

% At the moment we only implement the Siddon ray-tracer in [1]
hTracer = matRad_RayTracerSiddon();

if nargin >= 5
    hTracer.lateralCutOff = lateralCutoff;
end

if nargin >= 4
    [radDepthV,radDepthCube] = hTracer.traceCube(stfElement,ct,V,rot_coordsV);
elseif nargin >= 3
    [radDepthV,radDepthCube] = hTracer.traceCube(stfElement,ct,V);
else
    [radDepthV,radDepthCube] = hTracer.traceCube(stfElement,ct);
end

%{
    
%If no subset of voxels is specified, take all of them
if nargin < 3
    V = transpose(1:prod(ct.cubeDim));
end

%if we don't provide rotated patient coordinates we can compute
%them here on our own
if nargin < 4
    [yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);
    xCoordsV       = xCoordsV_vox(:)*ct.resolution.x-stfElement.isoCenter(1);
    yCoordsV       = yCoordsV_vox(:)*ct.resolution.y-stfElement.isoCenter(2);
    zCoordsV       = zCoordsV_vox(:)*ct.resolution.z-stfElement.isoCenter(3);
    coordsV        = [xCoordsV yCoordsV zCoordsV];
    
    % Get Rotation Matrix
    % Do not transpose matrix since we usage of row vectors &
    % transformation of the coordinate system need double transpose
    
    rotMat_system_T = matRad_getRotationMatrix(stfElement.gantryAngle,stfElement.couchAngle);
    
    % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
    rot_coordsV         = coordsV*rotMat_system_T;
    
    %
    rot_coordsV(:,1) = rot_coordsV(:,1)-stfElement.sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-stfElement.sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-stfElement.sourcePoint_bev(3);
end

%If we don't provide a lateral cutoff
if nargin < 5
    matRad_cfg = MatRad_Config.instance();
    lateralCutoff = matRad_cfg.propDoseCalc.defaultGeometricCutOff;
end

% set up rad depth cube for results
radDepthCube = repmat({NaN*ones(ct.cubeDim)},ct.numOfCtScen);

% set up ray matrix direct behind last voxel
rayMx_bev_y = max(rot_coordsV(:,2)) + max([ct.resolution.x ct.resolution.y ct.resolution.z]);
rayMx_bev_y = rayMx_bev_y + stfElement.sourcePoint_bev(2);

% set up list with bev coordinates for calculation of radiological depth
coords = zeros(prod(ct.cubeDim),3);
coords(V,:) = rot_coordsV;

% calculate spacing of rays on ray matrix
rayMxSpacing = 1/sqrt(2) * min([ct.resolution.x ct.resolution.y ct.resolution.z]);

% define candidate ray matrix covering 1000x1000mm^2
numOfCandidateRays = 2 * ceil(500/rayMxSpacing) + 1;
candidateRayMx     = zeros(numOfCandidateRays);

% define coordinates
[candidateRaysCoords_X,candidateRaysCoords_Z] = meshgrid(rayMxSpacing*[floor(-500/rayMxSpacing):ceil(500/rayMxSpacing)]);

% check which rays should be used
for i = 1:stfElement.numOfRays
    
    ix = (candidateRaysCoords_X(:) - (1+rayMx_bev_y/stfElement.SAD)*stfElement.ray(i).rayPos_bev(1)).^2 + ...
        (candidateRaysCoords_Z(:) - (1+rayMx_bev_y/stfElement.SAD)*stfElement.ray(i).rayPos_bev(3)).^2 ...
        <= lateralCutoff^2;
    
    candidateRayMx(ix) = 1;
    
end

% set up ray matrix
rayMx_bev = [candidateRaysCoords_X(logical(candidateRayMx(:))) ...
    rayMx_bev_y*ones(sum(candidateRayMx(:)),1) ...
    candidateRaysCoords_Z(logical(candidateRayMx(:)))];

%     figure,
%     for jj = 1:length(rayMx_bev)
%        plot(rayMx_bev(jj,1),rayMx_bev(jj,3),'rx'),hold on
%     end

% Rotation matrix. Transposed because of row vectors
rotMat_vectors_T = transpose(matRad_getRotationMatrix(stfElement.gantryAngle,stfElement.couchAngle));

% rotate ray matrix from bev to world coordinates
rayMx_world = rayMx_bev * rotMat_vectors_T;

% criterium for ray selection
raySelection = rayMxSpacing/2;

% perform ray tracing over all rays
for i = 1:size(rayMx_world,1)
    
    % run siddon ray tracing algorithm
    [~,l,rho,~,ixHitVoxel] = matRad_siddonRayTracer(stfElement.isoCenter, ...
        ct.resolution, ...
        stfElement.sourcePoint, ...
        rayMx_world(i,:), ...
        ct.cube);
    
    % find voxels for which we should remember this tracing because this is
    % the closest ray by projecting the voxel coordinates to the
    % intersection points with the ray matrix and checking if the distance
    % in x and z direction is smaller than the resolution of the ray matrix
    scale_factor = (rayMx_bev_y - stfElement.sourcePoint_bev(2)) ./ ...
        coords(ixHitVoxel,2);
    
    x_dist = coords(ixHitVoxel,1).*scale_factor - rayMx_bev(i,1);
    z_dist = coords(ixHitVoxel,3).*scale_factor - rayMx_bev(i,3);
    
    ixRememberFromCurrTracing = x_dist > -raySelection & x_dist <= raySelection ...
        & z_dist > -raySelection & z_dist <= raySelection;
    
    if any(ixRememberFromCurrTracing)
        
        for j = 1:ct.numOfCtScen
            % calc radiological depths
            
            % eq 14
            % It multiply voxel intersections with \rho values.
            d = l .* rho{j}; %Note. It is not a number "one"; it is the letter "l"
            
            % Calculate accumulated d sum.
            dCum = cumsum(d)-d/2;
            
            % write radiological depth for voxel which we want to remember
            radDepthCube{j}(ixHitVoxel(ixRememberFromCurrTracing)) = dCum(ixRememberFromCurrTracing);
        end
    end
    
end

% only take voxel inside the patient
for i = 1:ct.numOfCtScen
    radDepthV{i} = radDepthCube{i}(V);
end

%}

