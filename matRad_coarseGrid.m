function [Vcoarse,cubeDimCoarse,vXcoarse,vYcoarse,vZcoarse] = matRad_coarseGrid(ct,newResolution,V)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad generating a coarse/dense voxel grid
% 
% call
%   [Vcoarse,cubeDimCoarse,vXcoarse,vYcoarse,vZcoarse] = matRad_coarseGrid(ct,newResolution,V)
% input
%   ct:                    matRad ct structure
%   newResolution:         [x y z] (mm)
%   V:                     linear voxels indices from cst at original CT resolution
%
% output
%   Vcoarse:        voxels indexes from cst for coarse grid
%   cubeDimCoarse:  coarse cube dimensions
%   vXcoarse,vYcoarse,vZcoarse: 
%
% References
%   -
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

vXcoarse = ct.x(1):newResolution.x:ct.x(end);
vYcoarse = ct.y(1):newResolution.y:ct.y(end);
vZcoarse = ct.z(1):newResolution.z:ct.z(end);

tmpCube    = zeros(ct.cubeDim);
tmpCube(V) = 1;

[ Y,  X,  Z] = meshgrid(ct.x,ct.y,ct.z);
[Yq, Xq, Zq] = meshgrid(vXcoarse,vYcoarse,vZcoarse);

% interpolate cube - cube is now stored in Y X Z 
Vcoarse = find(interp3(Y,X,Z,tmpCube,Yq,Xq,Zq)>0.5);

cubeDimCoarse = [numel(vXcoarse) numel(vXcoarse) numel(vXcoarse)];


end

