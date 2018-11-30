function [radDepthVcoarse,radDepthIxcoarse] = matRad_interpRadDepth(ct,V,radDepthIx,radDepthV,vXgrid,vYgrid,vZgrid,Vcoarse)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% down/up sampling the radiological depth dose cubes
% 
% call
%   [radDepthVcoarse,radDepthIxcoarse] = matRad_interpRadDepth(ct,V,radDepthIx,radDepthV,vXgrid,vYgrid,vZgrid,Vcoarse)
% input
%   ct: matRad ct structure
%   V: voxel numbers 
%   radDepthIx: index in V of voxels with radiological depths 
%   radDepthV: V of cst voxels with radiological depths 
%   vXgrid,vYgrid,vZgrid: downsampled cube grid 
%   Vcoarse: downsampled V 
%
% output
%   radDepthVcoarse: V of cst voxels with radiological depths for coarse cube
%   radDepthIxcoarse: index in V of voxels with radiological depths for
%   coarse cube
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
% future: loop over multiple ct
radDepthCube                = NaN*ones(ct.cubeDim);
radDepthCube(V(radDepthIx)) = radDepthV{1}(radDepthIx);

[ Y,  X,  Z] = meshgrid(ct.x,ct.y,ct.z);
[Yq, Xq, Zq] = meshgrid(vXgrid,vYgrid,vZgrid);

% interpolate cube - cube is now stored in Y X Z 
coarseRadDepthCube   = (interp3(Y,X,Z,radDepthCube,Yq,Xq,Zq));
radDepthVcoarse{1}     = coarseRadDepthCube(Vcoarse);
radDepthlinIxcoarse  = find(~isnan(coarseRadDepthCube));
[~,radDepthIxcoarse] = ismember(radDepthlinIxcoarse,Vcoarse);
 
end

