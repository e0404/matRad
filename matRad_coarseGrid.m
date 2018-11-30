function [Vcoarse,cubeDimCoarse,vXcoarse,vYcoarse,vZcoarse] = matRad_coarseGrid(ct,NewResolution,V)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad generating a coarse/dense voxel grid
% 
% call
%   [Vcoarse,cubeDimCoarse,vXcoarse,vYcoarse,vZcoarse] = matRad_coarseGrid(ct,NewResolution,V)
% input
%   ct: matRad ct structure
%   NewResolution : [x y z] (mm)
%   V: voxels indexes from cst
%
% output
%   Vcoarse: voxels indexes from cst for coarse grid
%   cubeDimCoarse : coarse cube dimensions
%   vXcoarse,vYcoarse,vZcoarse: 
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vXcoarse = ct.x(1):NewResolution(1):ct.x(end);
vYcoarse = ct.y(1):NewResolution(2):ct.y(end);
vZcoarse = ct.z(1):NewResolution(3):ct.z(end);

tmpCube    = zeros(ct.cubeDim);
tmpCube(V) = 1;

[ Y,  X,  Z] = meshgrid(ct.x,ct.y,ct.z);
[Yq, Xq, Zq] = meshgrid(vXcoarse,vYcoarse,vZcoarse);

% interpolate cube - cube is now stored in Y X Z 
Vcoarse = find(interp3(Y,X,Z,tmpCube,Yq,Xq,Zq));

cubeDimCoarse = [numel(vXcoarse) numel(vXcoarse) numel(vXcoarse)];


end

