function [stf, dij] = matRad_addMUdataFromMachine(machine, stf, dij)
% helper function to add MU data included in machine file to dij and stf
% computed with earlier versions of matRad
% 
% call
%   [stf, dij] = matRad_addMUdataFromMachine(machine, stf, dij)
%   stf = matRad_addMUdataFromMachine(machine, stf)
%
% input
%   machine:        machine struct as stored in basedata file
%   stf:            matRad steering information struct
%   dij:            matRad dose influence matrix struct
%
% output
%   stf:            matRad steering information struct with MUdata
%   dij:            matRad dose influence matrix struct with MUdata
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2012 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Each base data entry should have a field called MUdata which contains
%minMU, maxMU and numParticlesPerMU
MUdata             = [machine.data.MUdata];

bdMinMU = [MUdata.minMU];
bdMaxMU = [MUdata.maxMU];
bdNumParticlesPerMU = [MUdata.numParticlesPerMU];
bdEnergies    = [machine.data.energy];

%Adding info to stf file
for iBeam = 1:numel(stf)
    for iRay = 1:stf(iBeam).numOfRays
        energy         = stf(iBeam).ray(iRay).energy;
        minMUSpot = spline(bdEnergies,bdMinMU,energy); %Will change this
        maxMUSpot = spline(bdEnergies,bdMaxMU,energy);
        NumParticlesPerMUspot      = spline(bdEnergies,bdNumParticlesPerMU,energy);

        [stf.ray(iRay).minMU] = minMUSpot;
        [stf.ray(iRay).maxMU] = maxMUSpot;
        [stf.ray(iRay).numParticlesPerMU]   = NumParticlesPerMUspot;

    end
end

if nargin > 2 && nargout == 2
    %Adding info to dij
    dij.minMU               = zeros(dij.totalNumOfBixels,1);
    dij.maxMU               = Inf(dij.totalNumOfBixels,1);
    dij.numParticlesPerMU   = 1e6*ones(dij.totalNumOfBixels,1);
    for iBixel = 1:dij.totalNumOfBixels
        
        beamNum = dij.beamNum(iBixel);
        rayNum = dij.rayNum(iBixel);
        spotNum = dij.bixelNum(iBixel);

        minMUSpot   = stf(beamNum).ray(rayNum).minMU(spotNum);
        maxMUSpot   = stf(beamNum).ray(rayNum).maxMU(spotNum);
        numParticlesPerMUspot = stf(beamNum).ray(rayNum).numParticlesPerMU(spotNum);

        dij.minMU(iBixel,1) = minMUSpot;
        dij.mxaMU(iBixel,1) = maxMUSpot;
        dij.numParticlesPerMU(iBixel)   = numParticlesPerMUspot;
    end
end







