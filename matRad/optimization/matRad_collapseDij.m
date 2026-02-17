function dijNew = matRad_collapseDij(dij, mode)
% matRad collapse dij function for simulation of 3D conformal treatments.
% Function to supress intensity-modulation for photons in order to simulate 
% 3D conformal treatments.
%
% call
%   dijNew = matRad_collapseDij(dij)
%
% input
%   dij:    dose influence matrix
%   mode:   collpase mode, beam or ray
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

 if nargin < 2
     mode = 'beam'; % default
 end
 validatestring(mode, {'beam','ray'});

  switch mode
      case 'beam'
            dijNew.totalNumOfBixels = dij.numOfBeams;
            dijNew.totalNumOfRays   = dij.numOfBeams;
            dijNew.numOfBeams       = dij.numOfBeams;
            dijNew.numOfRaysPerBeam = ones(dij.numOfBeams,1);
            
            dijNew.beamNum          = (1:dij.numOfBeams)';
            dijNew.bixelNum         = ones(dij.numOfBeams, 1);
            dijNew.rayNum           = ones(dij.numOfBeams,1);
      case 'ray'
            dijNew.totalNumOfBixels = dij.totalNumOfRays;
            dijNew.totalNumOfRays   = dij.totalNumOfRays;
            dijNew.numOfBeams       = dij.numOfBeams;
            dijNew.numOfRaysPerBeam = dij.numOfRaysPerBeam;

            dijNew.beamNum          = repelem(1:dij.numOfBeams, dij.numOfRaysPerBeam(:))';
            dijNew.bixelNum         = ones(dij.totalNumOfRays, 1);
            dijNew.rayNum           = cell2mat(arrayfun(@(n) 1:n, dij.numOfRaysPerBeam, 'UniformOutput', false))';

            totNumRays =  [0,cumsum(dij.numOfRaysPerBeam)];
  end

if isfield(dij,'numParticlesPerMU')
    switch mode
        case 'beam'
            dijNew.numParticlesPerMU = zeros(dij.numOfBeams,1);
            for j = 1:dij.numOfBeams
                dijNew.numParticlesPerMU(j) = sum(dij.numParticlesPerMU(dij.beamNum == j));
            end
        case 'ray'
            dijNew.numParticlesPerMU = zeros(dij.totalNumOfRays,1);
            for j = 1:dij.numOfBeams
                for k = 1:dij.numOfRaysPerBeam(j)
                   dijNew.numParticlesPerMU(totNumRays(j) + k) = sum(dij.numParticlesPerMU((dij.rayNum == k)&(dij.beamNum == j)));
                end
            end
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
        switch mode
            case 'beam'
                tmp = sparse(dij.doseGrid.numOfVoxels,dij.numOfBeams);          % initialize sparse matrix
                for j = 1:dij.numOfBeams
                    % Sum only the columns corresponding to beam j
                    tmp(:, j) = sum(dij.(quantityName){i}(:, dij.beamNum == j), 2);
                end  
            case 'ray'
                tmp = sparse(dij.doseGrid.numOfVoxels,dij.totalNumOfRays);
                for j = 1:dij.numOfBeams
                    for k = 1:dij.numOfRaysPerBeam(j)
                        tmp(:, totNumRays(j) + k) = sum(dij.(quantityName){i}(:, (dij.rayNum == k)&(dij.beamNum == j)),2);
                    end
                end
        end
        dijNew.(quantityName){i} = tmp;
     end

end

