function voxelProbCube = matRad_sampleVoxelProb(cst,ct,sigma,VOIName,ncase)

rng(0);

% get cst index if VOI
cstidx = find(strcmp([cst(:,2)],VOIName));

% create empty Cube
voxelProbCube = zeros(ct.cubeDim);

% create samples
xShift_vox = round(sigma(1).*randn(1,ncase)./ct.resolution.x);
yShift_vox = round(sigma(2).*randn(1,ncase)./ct.resolution.y);
zShift_vox = round(sigma(3).*randn(1,ncase)./ct.resolution.z);

% get voxel cooridinates
[yCoordsVOI_vox, xCoordsVOI_vox, zCoordsVOI_vox] = ind2sub(ct.cubeDim,cst{cstidx,4}{1});

% loop over all samples
for i = 1:ncase
            
    shiftedVOIidx = sub2ind(ct.cubeDim, yCoordsVOI_vox + yShift_vox(i),...
                                        xCoordsVOI_vox + xShift_vox(i),...
                                        zCoordsVOI_vox + zShift_vox(i));
                                     
                                     
    voxelProbCube(shiftedVOIidx) = voxelProbCube(shiftedVOIidx) + 1/ncase;
    
end

end

