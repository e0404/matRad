function [minAbsDist,minDist] = matRad_calcMinDist(ct,VOI1VocelIDs,VOI2VoxelIDs)

    % get cube coordinates
    [xCoords_vox, yCoords_vox, zCoords_vox] = meshgrid(1:ct.cubeDim(2),1:ct.cubeDim(1),1:ct.cubeDim(3));

    xCoords = xCoords_vox(:)*ct.resolution.x;
    yCoords = yCoords_vox(:)*ct.resolution.y;
    zCoords = zCoords_vox(:)*ct.resolution.z;

    % get surface voxel of VOI2
    VOI2SurfaceVoxelIDs = matRad_getVOISurfaceVoxel(ct,VOI2VoxelIDs);

    % get VOI coordinates (consider only surface voxels for VOI2)
    VOI1coords = [xCoords(VOI1VocelIDs) yCoords(VOI1VocelIDs) zCoords(VOI1VocelIDs)];
    VOI2coords = [xCoords(VOI2SurfaceVoxelIDs) yCoords(VOI2SurfaceVoxelIDs) zCoords(VOI2SurfaceVoxelIDs)];

    % calculate min distance of VOI1 voxels to VOI2 voxels
    for i = 1:size(VOI1coords,1)
        
        % calculate absolute min distance
        coordsDifference    = repmat(VOI1coords(i,:),size(VOI2coords,1),1) - VOI2coords;
        [minAbsDist(i),idx] = min(sqrt(sum((coordsDifference).^2,2)));
        
        % divide absolute min distance in x,y and z
        minDist(:,i)        = abs(coordsDifference(idx,:));
        
    end

end
