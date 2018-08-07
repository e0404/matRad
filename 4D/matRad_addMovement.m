function ct = matRad_addMovement(ct, motionPeriod, numOfCtScen, amp)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adds artificial patient movement by creating a Deformation Vector Filed
% and applying it to the ct.cube by geometric transformation of the ct
%
% call
%   ct = matRad_addmovement(ct, ct.motionPeriod, ct.numOfCtScen, amp)
%
% input
%   ct :            ct cube
%   motionPeriod:   the extent of a whole breathing cycle (in seconds)
%   numOfCtScen:    number of the desired phases
%   amp:            amplitude of the sinosoidal movement(in pixels)
%
% output
%   ct :            dvf and cube fields for each scenario added
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

if nargin < 4
    amp = [0 5 0];
end

% if there is already a dvf, no need to addMovement
if isfield(ct,'dvf') 
    return
end

ct.motionPeriod = motionPeriod;
ct.numOfCtScen = numOfCtScen;

for i = 1:numOfCtScen
    
    im = ct.cube{1};
    
    ct.dvf{i} = zeros([size(im), 3]);
    
    ct.dvf{i}(:,:,:,1) = amp(1) * sin((i - 1) * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,2) = amp(2) * sin((i - 1) * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,3) = amp(3) * sin((i - 1) * pi / (motionPeriod * 2));
    

    ct.cube{i} = imwarp(im, ct.dvf{i});
    ct.cubeHU{i} = imwarp(im, ct.dvf{i});
    
    ct.dvf{i}(:,:,:,1) = ct.dvf{i}(:,:,:,1) * ct.resolution.x;
    ct.dvf{i}(:,:,:,2) = ct.dvf{i}(:,:,:,2) * ct.resolution.y;
    ct.dvf{i}(:,:,:,3) = ct.dvf{i}(:,:,:,3) * ct.resolution.z;
    ct.dvf{i} = permute(ct.dvf{i}, [4,1,2,3]);
    
end
end

