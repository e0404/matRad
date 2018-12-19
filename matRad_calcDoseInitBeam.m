fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);

% remember beam and bixel number
if calcDoseDirect
    dij.beamNum(i)    = i;
    dij.rayNum(i)     = i;
    dij.bixelNum(i)   = i;
end

bixelsPerBeam = 0;

% convert voxel indices to real coordinates using iso center of beam i
xCoordsV       = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
yCoordsV       = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
zCoordsV       = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
coordsV        = [xCoordsV yCoordsV zCoordsV];

xCoordsVdoseGrid = xCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.x-stf(i).isoCenter(1);
yCoordsVdoseGrid = yCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.y-stf(i).isoCenter(2);
zCoordsVdoseGrid = zCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.z-stf(i).isoCenter(3);
coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];

% Get Rotation Matrix
% Do not transpose matrix since we usage of row vectors &
% transformation of the coordinate system need double transpose

rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

% Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
rot_coordsV         = coordsV*rotMat_system_T;
rot_coordsVdoseGrid = coordsVdoseGrid*rotMat_system_T;

rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

rot_coordsVdoseGrid(:,1) = rot_coordsVdoseGrid(:,1)-stf(i).sourcePoint_bev(1);
rot_coordsVdoseGrid(:,2) = rot_coordsVdoseGrid(:,2)-stf(i).sourcePoint_bev(2);
rot_coordsVdoseGrid(:,3) = rot_coordsVdoseGrid(:,3)-stf(i).sourcePoint_bev(3);

% calculate geometric distances
geoDistVdoseGrid{1}= sqrt(sum(rot_coordsVdoseGrid.^2,2));

% Calculate radiological depth cube
fprintf('matRad: calculate radiological depth cube...');
radDepthVctGrid = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,effectiveLateralCutoff);
fprintf('done.\n');

% interpolate radiological depth cube to dose grid resolution
radDepthVdoseGrid = matRad_interpRadDepth...
    (ct,1,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);
 
% limit rotated coordinates to positions where ray tracing is availabe
rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(radDepthVdoseGrid{1}),:);
 