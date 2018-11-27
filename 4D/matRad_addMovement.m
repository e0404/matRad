function [ct, cst] = matRad_addMovement(ct, cst, motionPeriod, numOfCtScen, amp,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adds artificial sinosodal patient motion by creating a deformation vector
% field and applying it to the ct.cube by geometric transformation
%
% call
%   ct = matRad_addMovement(ct, ct.motionPeriod, ct.numOfCtScen, amp)
%
% input
%   ct:             matRad ct struct
%   cst:            matRad cst struct
%   motionPeriod:   the length of a whole breathing cycle (in seconds)
%   numOfCtScen:    number of ct phases
%   amp:            amplitude of the sinosoidal movement (in pixels)
%   visBool         boolean flag for visualization
%
%   note:           1st dim --> x LPS coordinate system
%                   2nd dim --> y LPS coordinate system
%                   3rd dim --> z LPS coordinate system
%                   a positive amplitude moves the phantom to the right,
%                   anterior, inferior
%
% output
%   ct:             modified matRad ct struct including dvf and cubes for 
%                   all phases
%   cst:            modified matRad cst struct
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
    
    ct.dvf{i}(:,:,:,1) = amp(1) * sin((i-1)*pi / numOfCtScen)^2; % deformation along x direction (i.e. 2nd coordinate in dose/ct)
    ct.dvf{i}(:,:,:,2) = amp(2) * sin((i-1)*pi / numOfCtScen)^2; % deformation along y direction (i.e. 3rd coordiatie in dose/ct)
    ct.dvf{i}(:,:,:,3) = amp(3) * sin((i-1)*pi / numOfCtScen)^2;
    
    % warp ct
    if isfield(ct,'hlut')
       padValue = min(ct.hlut(:,2));
    else
       padValue = -1024;
    end
    
    ct.cubeHU{i} = imwarp(ct.cubeHU{1}, ct.dvf{i},'FillValues',padValue);
    
    if isfield(ct,'cube')
       ct.cube{i}   = imwarp(ct.cube{1},   ct.dvf{i},'FillValues',0);
    end
    
 
    
    % warp cst
    for j = 1:size(cst,1)
        tmp = zeros(ct.cubeDim);
        tmp(cst{j,4}{1}) = 1;
        tmpWarp     = imwarp(tmp, ct.dvf{i});
        
        cst{j,4}{i} = find(tmpWarp > .5);
    end
    
    % convert dvfs to [mm]
    tmp = ct.dvf{i}(:,:,:,1);
    ct.dvf{i}(:,:,:,1) = -ct.dvf{i}(:,:,:,2) * ct.resolution.x;
    ct.dvf{i}(:,:,:,2) = -tmp * ct.resolution.y;
    ct.dvf{i}(:,:,:,3) = -ct.dvf{i}(:,:,:,3) * ct.resolution.z;
    
    ct.dvf{i} = permute(ct.dvf{i}, [4,1,2,3]);
    
end



if visBool
   slice       = round(ct.cubeDim(3)/2);
   figure,
   for i = 1:numOfCtScen
      clf,
      imagesc(ct.cubeHU{i}(:,:,slice))
      pause(.5);
   end
end





