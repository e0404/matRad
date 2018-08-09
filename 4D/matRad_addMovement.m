function ct = matRad_addMovement(ct, motionPeriod, numOfCtScen, amp)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adds artificial sinosodal patient motion by creating a deformation vector
% field and applying it to the ct.cube by geometric transformation
%
% call
%   ct = matRad_addMovement(ct, ct.motionPeriod, ct.numOfCtScen, amp)
%
% input
%   ct:             matRad ct struct
%   motionPeriod:   the length of a whole breathing cycle (in seconds)
%   numOfCtScen:    number of ct phases
%   amp:            amplitude of the sinosoidal movement (in pixels)
%
% output
%   ct :            modified matRad ct struct including dvf and cubes for 
%                   all phases
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% book keeping
ct.motionPeriod = motionPeriod;
ct.numOfCtScen = numOfCtScen;

% set type
ct.dvfType = 'pull'; % push or pull

% generate scenarios
for i = 1:numOfCtScen
        
    ct.dvf{i} = zeros([ct.cubeDim, 3]);
    
    ct.dvf{i}(:,:,:,1) = amp(1) * sin((i-1)*pi / numOfCtScen)^2;
    ct.dvf{i}(:,:,:,2) = amp(2) * sin((i-1)*pi / numOfCtScen)^2;
    ct.dvf{i}(:,:,:,3) = amp(3) * sin((i-1)*pi / numOfCtScen)^2;
    
    ct.cube{i}   = imwarp(ct.cube{1},   ct.dvf{i});
    ct.cubeHU{i} = imwarp(ct.cubeHU{1}, ct.dvf{i});
    
    ct.dvf{i}(:,:,:,1) = ct.dvf{i}(:,:,:,1) * ct.resolution.x;
    ct.dvf{i}(:,:,:,2) = ct.dvf{i}(:,:,:,2) * ct.resolution.y;
    ct.dvf{i}(:,:,:,3) = ct.dvf{i}(:,:,:,3) * ct.resolution.z;
    
    ct.dvf{i} = permute(ct.dvf{i}, [4,1,2,3]);
    
end
