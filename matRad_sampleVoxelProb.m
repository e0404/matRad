function [voxelProbCube,voxelShift] = matRad_sampleVoxelProb(cst,ct,sigma,VOIName,ncase)

rng(0);

% get cst index if VOI
cstidx = find(strcmp([cst(:,2)],VOIName));

% create empty Cube
voxelProbCube = zeros(ct.cubeDim);

% create samples
voxelXShift = round(sigma(1).*randn(1,ncase)./ct.resolution.x);
voxelYShift = round(sigma(2).*randn(1,ncase)./ct.resolution.y);
voxelZShift = round(sigma(3).*randn(1,ncase)./ct.resolution.z);
voxelShift  = [voxelXShift;voxelYShift;voxelZShift];

% get voxel cooridinates
[yCoordsVOI_vox, xCoordsVOI_vox, zCoordsVOI_vox] = ind2sub(ct.cubeDim,cst{cstidx,4}{1});

% loop over all samples
for i = 1:ncase
            
    shiftedVOIidx = sub2ind(ct.cubeDim, yCoordsVOI_vox + voxelYShift(i),...
                                        xCoordsVOI_vox + voxelXShift(i),...
                                        zCoordsVOI_vox + voxelZShift(i));
                                                                       
    voxelProbCube(shiftedVOIidx) = voxelProbCube(shiftedVOIidx) + 1/ncase;
    
end

end

