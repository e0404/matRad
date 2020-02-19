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
% [radDepthCubeCtGrid, radDepthVctGrid] = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,effectiveLateralCutoff);

[~, energyIx] =  intersect([machine.data(:).energy], stf.ray.energy);

yd = diff(machine.data(energyIx).initFocus.sigma)./diff(machine.data(energyIx).initFocus.dist);
xd = (machine.data(energyIx).initFocus.dist(2:end)+machine.data(energyIx).initFocus.dist(1:(end-1)))/2;

initSigmaConv = interp1(xd,yd,stf.ray.SSD);

if strcmp(anaMode, 'stdCorr')
    kernelSize = min(ct.cubeDim(2), ct.cubeDim(3)) + (mod(min(ct.cubeDim(2), ct.cubeDim(3)),2) -  1);
    res = [ct.resolution.y, ct.resolution.z];
    dist = 0;
    for ii = 1:ct.cubeDim(1)
        slice = reshape(radDepthCubeCtGrid(ii,:,:),ct.cubeDim(2),ct.cubeDim(3));
        notNaN = find(~isnan(slice));
        slice(isnan(slice)) = 0;
        dist = sum(slice(notNaN))/ numel(slice(notNaN));
        
        std1(ii) = std(slice,0,'all')^2 + 1;
        dist11(ii) = sum(1./(std1)) / 2 * pi;
        
        nPadd = floor(kernelSize/2);
        paddedSlice = zeros(size(slice,1) + 2 * nPadd, size(slice,2) + 2 * nPadd);
        paddedSlice(nPadd + 1:ct.cubeDim(2) + nPadd, nPadd + 1:ct.cubeDim(3) + nPadd) = slice;
        
        paddedSlice(1:nPadd, 1:nPadd) = slice(1,1);
        paddedSlice(ct.cubeDim(2) + nPadd + 1:end, 1:nPadd) = slice(end,1);
        paddedSlice(1:nPadd, ct.cubeDim(3) + nPadd + 1:end) = slice(1,end);
        paddedSlice(ct.cubeDim(2) + nPadd + 1:end, ct.cubeDim(3) + nPadd + 1:end) = slice(end);
        
        paddedSlice(1:nPadd,nPadd + 1:ct.cubeDim(3) + nPadd) = repmat(slice(1,:),nPadd,1);
        paddedSlice(ct.cubeDim(2) + nPadd + 1:end,nPadd + 1:ct.cubeDim(3) + nPadd) = repmat(slice(end,:),nPadd,1);
        paddedSlice(nPadd + 1:ct.cubeDim(2) + nPadd,1:nPadd) = repmat(slice(:,1),1,nPadd);
        paddedSlice(nPadd + 1:ct.cubeDim(2) + nPadd ,ct.cubeDim(3) + nPadd + 1:end) = repmat(slice(:,end),1,nPadd);



        if isnan(dist)
            sigma = 0;
        else
%             sigma = testTMP/150 * dist;
            sigma = dist11(ii);
        end
        kernel = matRad_create2dimGaussKernel(kernelSize, sigma, res);

        tmpCstd     = (convn(paddedSlice.^2, kernel, 'same') - convn(paddedSlice, kernel, 'same').^2);
        tmpCstd(tmpCstd < 0) = 0;
        tmpCstd = sqrt(tmpCstd);
        tmpRadDepth = convn(paddedSlice, kernel, 'same');
%         imagesc(tmpCstd);
%         figure
%         imagesc(tmpRadDepth);
%         figure
%         imagesc(slice)
%         close all
        cStdCtGrid(ii,:,:)   = tmpCstd((kernelSize+1)/2:size(tmpCstd,1)-(kernelSize-1)/2,(kernelSize+1)/2:size(tmpCstd,2)-(kernelSize-1)/2);
        meanRadDepths(ii,:,:) = tmpRadDepth((kernelSize+1)/2:size(tmpRadDepth,1)-(kernelSize-1)/2,(kernelSize+1)/2:size(tmpRadDepth,2)-(kernelSize-1)/2);    
    end

%     assignin('base','cStdCtGrid',cStdCtGrid);
%     assignin('base','meanRadDepths',meanRadDepths);
%     assignin('base','radDepthCubeCtGrid',radDepthCubeCtGrid);
    
    cStdCtGrid(isnan(cStdCtGrid)) = 0;
    cStdVctGrid = {cStdCtGrid(VctGrid)};
    
  
    meanRadDepths = meanRadDepths(VctGrid);
    
    cSstVdoseGrid = matRad_interpRadDepth...
    (ct,1,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,cStdVctGrid);
end

fprintf('done.\n');

% interpolate radiological depth cube and std cube to dose grid resolution
radDepthVdoseGrid = matRad_interpRadDepth...
    (ct,1,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);

% limit rotated coordinates to positions where ray tracing is availabe
rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(radDepthVdoseGrid{1}),:);
 