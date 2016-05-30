function [voxelProbCube,voxelShift,idxShift] = matRad_sampleVoxelProb(cst,ct,sigma,VOIName,ncase)

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

% convert voxel shifts to index shifts
idxShift    = voxelShift(2,:) + voxelShift(1,:).*ct.cubeDim(1) + voxelShift(3,:).*ct.cubeDim(2).*ct.cubeDim(1);

% loop over all samples
for i = 1:ncase
    
    shiftedVOIidx                = cst{cstidx,4}{1} + idxShift(i);
    voxelProbCube(shiftedVOIidx) = voxelProbCube(shiftedVOIidx) + 1/ncase;
    
end

end

