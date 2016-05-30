function [voxelProbCube,voxelShift,idxShift] = matRad_sampleVoxelProb(cst,ct,sigma,VOIName,ncase)

rng(0);

% get cst index if VOI
cstidx = find(strcmp([cst(:,2)],VOIName));

% create empty Cube
voxelProbCube = zeros(ct.cubeDim);

% create samples
voxelXShift = [0,round(sigma(1).*randn(1,ncase-1)./ct.resolution.x)];
voxelYShift = [0,round(sigma(2).*randn(1,ncase-1)./ct.resolution.y)];
voxelZShift = [0,round(sigma(3).*randn(1,ncase-1)./ct.resolution.z)];
voxelShift  = [voxelXShift;voxelYShift;voxelZShift];

% convert voxel shifts to index shifts
idxShift    = voxelShift(2,:) + voxelShift(1,:).*ct.cubeDim(1) + voxelShift(3,:).*ct.cubeDim(2).*ct.cubeDim(1);

% loop over all samples
for i = 1:ncase
    
    shiftedVOIidx                = cst{cstidx,4}{1} + idxShift(i);
    voxelProbCube(shiftedVOIidx) = voxelProbCube(shiftedVOIidx) + 1/ncase;
    
end

end

