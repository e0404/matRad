function dijNew = matRad_collapseDij(dij)
% matRad collapse dij function to to supress intensity-modulation for 
% photons in order to simulate 3D conformal treatments.
%
% call
%   dijNew = matRad_collapseDij(dij)
%
% input
%   dij:            dose influence matrix
%
% output
%   dij:            collapsed dose influence matrix
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dijNew.totalNumOfBixels = 1;
dijNew.totalNumOfRays   = 1;
dijNew.numOfBeams       = 1;
dijNew.numOfRaysPerBeam = 1;

dijNew.beamNum          = 1;
dijNew.bixelNum         = 1;
dijNew.rayNum           = 1;

dijNew.doseGrid = dij.doseGrid;
dijNew.ctGrid   = dij.ctGrid;

dijNew.numOfScenarios = dij.numOfScenarios;

for i = 1:dij.numOfScenarios
    dijNew.physicalDose{i} = sum(dij.physicalDose{i},2);
end

