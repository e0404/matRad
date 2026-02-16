function dijNew = matRad_collapseDij(dij)
% matRad collapse dij function for simulation of 3D conformal treatments.
% Function to supress intensity-modulation for photons in order to simulate 
% 3D conformal treatments.
%
% call
%   dijNew = matRad_collapseDij(dij)
%
% input
%   dij:    dose influence matrix
%
% output
%   dijNew: collapsed dose influence matrix
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dijNew.totalNumOfBixels = dij.numOfBeams;
dijNew.totalNumOfRays   = dij.numOfBeams;
dijNew.numOfBeams       = dij.numOfBeams;
dijNew.numOfRaysPerBeam = ones(dij.numOfBeams,1);

dijNew.beamNum          = (1:dij.numOfBeams)';
dijNew.bixelNum         = ones(dij.numOfBeams, 1);
dijNew.rayNum           = ones(dij.numOfBeams,1);

if isfield(dij,'numParticlesPerMU')
    dijNew.numParticlesPerMU = zeros(dij.numOfBeams,1);
    for j = 1:dij.numOfBeams
        dijNew.numParticlesPerMU(j) = sum(dij.numParticlesPerMU(dij.beamNum == j));
    end
end

dijNew.doseGrid = dij.doseGrid;
dijNew.ctGrid   = dij.ctGrid;

dijNew.numOfScenarios = dij.numOfScenarios;

collapsableQuantites = {'physicalDose', 'mLETDose', 'mAlphaDose', 'mSqrtBetaDose'};

% Identify quantities present in dij
quantitiesToCollapse = collapsableQuantites(ismember(collapsableQuantites, fieldnames(dij)));

for q = 1:numel(quantitiesToCollapse)
    quantityName = quantitiesToCollapse{q};
    
    for i = 1:numel(dij.(quantityName))
        if isempty(dij.(quantityName){i})
            dijNew.(quantityName){i} = [];
            continue;
        end

        tmp = sparse(dij.doseGrid.numOfVoxels,dij.numOfBeams);          % initialize sparse matrix
        for j = 1:dij.numOfBeams
            % Sum only the columns corresponding to beam j
            tmp(:, j) = sum(dij.(quantityName){i}(:, dij.beamNum == j), 2);
        end     
        dijNew.(quantityName){i} = tmp;
    end

end

