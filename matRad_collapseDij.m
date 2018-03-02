function dijNew = matRad_collapseDij(dij)

dijNew.totalNumOfBixels = 1;
dijNew.totalNumOfRays   = 1;
dijNew.numOfBeams       = 1;
dijNew.numOfRaysPerBeam = 1;

dijNew.beamNum          = 1;
dijNew.bixelNum         = 1;
dijNew.rayNum           = 1;

dijNew.numOfVoxels    = dij.numOfVoxels;
dijNew.resolution     = dij.resolution;
dijNew.dimensions     = dij.dimensions;
dijNew.numOfScenarios = dij.numOfScenarios;

for i = 1:dij.numOfScenarios
    dijNew.physicalDose{i} = sum(dij.physicalDose{i},2);
end