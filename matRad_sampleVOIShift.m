function VOIShift = matRad_sampleVOIShift(cst,ct,sigma,VOIName,ncase)

rng(0);

% get cst index if VOI
cstidx = find(strcmp([cst(:,2)],VOIName));

% create empty Cube
voxelProbCube = zeros(ct.cubeDim);

% create exact shifts
voxelXShift                = [0,sigma(1).*randn(1,ncase-1)./ct.resolution.x];
voxelYShift                = [0,sigma(2).*randn(1,ncase-1)./ct.resolution.y];
voxelZShift                = [0,sigma(3).*randn(1,ncase-1)./ct.resolution.z];
linInterpShift.voxelShift  = [voxelXShift;voxelYShift;voxelZShift];

voxelShiftX0 = floor(voxelXShift);
voxelShiftX1 = ceil(voxelXShift);
voxelShiftY0 = floor(voxelYShift);
voxelShiftY1 = ceil(voxelYShift);
voxelShiftZ0 = floor(voxelZShift);
voxelShiftZ1 = ceil(voxelZShift);

linInterpShift.idxShift.X0Y0Z0 = voxelShiftY0 + voxelShiftX0.*ct.cubeDim(1) + voxelShiftZ0.*ct.cubeDim(2).*ct.cubeDim(1);
linInterpShift.idxShift.X0Y0Z1 = voxelShiftY0 + voxelShiftX0.*ct.cubeDim(1) + voxelShiftZ1.*ct.cubeDim(2).*ct.cubeDim(1);
linInterpShift.idxShift.X0Y1Z0 = voxelShiftY1 + voxelShiftX0.*ct.cubeDim(1) + voxelShiftZ0.*ct.cubeDim(2).*ct.cubeDim(1);
linInterpShift.idxShift.X0Y1Z1 = voxelShiftY1 + voxelShiftX0.*ct.cubeDim(1) + voxelShiftZ1.*ct.cubeDim(2).*ct.cubeDim(1);
linInterpShift.idxShift.X1Y0Z0 = voxelShiftY0 + voxelShiftX1.*ct.cubeDim(1) + voxelShiftZ0.*ct.cubeDim(2).*ct.cubeDim(1);
linInterpShift.idxShift.X1Y0Z1 = voxelShiftY0 + voxelShiftX1.*ct.cubeDim(1) + voxelShiftZ1.*ct.cubeDim(2).*ct.cubeDim(1);
linInterpShift.idxShift.X1Y1Z0 = voxelShiftY1 + voxelShiftX1.*ct.cubeDim(1) + voxelShiftZ0.*ct.cubeDim(2).*ct.cubeDim(1);
linInterpShift.idxShift.X1Y1Z1 = voxelShiftY1 + voxelShiftX1.*ct.cubeDim(1) + voxelShiftZ1.*ct.cubeDim(2).*ct.cubeDim(1);

linInterpShift.idxShift.x = (voxelXShift - voxelShiftX0)./(voxelShiftX1 - voxelShiftX0);
linInterpShift.idxShift.y = (voxelYShift - voxelShiftY0)./(voxelShiftY1 - voxelShiftY0);
linInterpShift.idxShift.z = (voxelZShift - voxelShiftZ0)./(voxelShiftZ1 - voxelShiftZ0);

linInterpShift.idxShift.x(isnan(linInterpShift.idxShift.x)) = 0;
linInterpShift.idxShift.y(isnan(linInterpShift.idxShift.y)) = 0;
linInterpShift.idxShift.z(isnan(linInterpShift.idxShift.z)) = 0;

% round shifts to voxel dimension
voxelXShiftRound        = round(voxelXShift);
voxelYShiftRound        = round(voxelYShift);
voxelZShiftRound        = round(voxelZShift);
roundedShift.voxelShift = [voxelXShiftRound;voxelYShiftRound;voxelZShiftRound];
roundedShift.idxShift   = voxelYShiftRound + voxelXShiftRound.*ct.cubeDim(1) + voxelZShiftRound.*ct.cubeDim(2).*ct.cubeDim(1);

% use rounded shifts to calculate probability cube
for i = 1:ncase
    
    shiftedVOIidx                = cst{cstidx,4}{1} - roundedShift.idxShift(i);
    voxelProbCube(shiftedVOIidx) = voxelProbCube(shiftedVOIidx) + 1/ncase;
    
end

VOIShift.voxelProbCube  = voxelProbCube;
VOIShift.linInterpShift = linInterpShift;
VOIShift.roundedShift   = roundedShift;
VOIShift.ncase          = ncase;
VOIShift.shiftType      = 'rounded'; % default option for optimization

end

