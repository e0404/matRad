function minDist = matRad_calcMinDist(cst,ct,VOIName1,VOIName2)

    % get cst index
    cstidxVOI1 = find(~cellfun('isempty',strfind(cst(:,2),VOIName1)));
    cstidxVOI2 = find(~cellfun('isempty',strfind(cst(:,2),VOIName2)));

    % get cube coordinates
    [xCoords_vox, yCoords_vox, zCoords_vox] = meshgrid(1:ct.cubeDim(1),1:ct.cubeDim(2),1:ct.cubeDim(3));

    xCoords = xCoords_vox(:)*ct.resolution.x;
    yCoords = yCoords_vox(:)*ct.resolution.y;
    zCoords = zCoords_vox(:)*ct.resolution.z;

    % get surface voxel of VOI2
    VOI2SurfaceVoxel = matRad_getVOISurfaceVoxel(cst,ct,VOIName2);

    % get VOI coordinates (consider only surface voxels for VOI2)
    VOI1coords = [xCoords(cst{cstidxVOI1,4}{1}) yCoords(cst{cstidxVOI1,4}{1}) zCoords(cst{cstidxVOI1,4}{1})];
    VOI2coords = [xCoords(VOI2SurfaceVoxel) yCoords(VOI2SurfaceVoxel) zCoords(VOI2SurfaceVoxel)];

    % calculate min distance of VOI1 voxels to VOI2 voxels
    for i = 1:size(VOI1coords,1)
        
        distancesToVOI2 = sqrt(sum((repmat(VOI1coords(i,:),size(VOI2coords,1),1) - VOI2coords).^2,2));
        minDist(i)      = min(distancesToVOI2);
        
    end

end