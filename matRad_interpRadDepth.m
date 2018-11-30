function radDepthVcoarse = matRad_interpRadDepth(ct,ctScenNum,V,Vcoarse,vXgrid,vYgrid,vZgrid,radDepthV)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% down/up sampling the radiological depth dose cubes
% 
% call
%   [radDepthVcoarse,radDepthIxcoarse] = matRad_interpRadDepth(ct,V,radDepthIx,radDepthV,vXgrid,vYgrid,vZgrid,Vcoarse)
% input
%   ct:             matRad ct structure
%   V:              linear voxel indices of the cst 
%   radDepthV:      radiological depth of radDepthIx
%   vXgrid:         query points of now location in x dimension
%   vYgrid:         query points of now location in y dimension
%   vZgrid:         query points of now location in z dimension
%   Vcoarse:        linear voxel indices of the down sampled grid resolution
%
% output
%   radDepthVcoarse:   interpolated radiological depth of radDepthIx
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

[ Y,  X,  Z] = meshgrid(ct.x,ct.y,ct.z);
[Yq, Xq, Zq] = meshgrid(vXgrid,vYgrid,vZgrid);

radDepthCube                = NaN*ones(ct.cubeDim);
radDepthCube(V(find(~isnan(radDepthV{1})))) = radDepthV{ctScenNum}(find(~isnan(radDepthV{1})));

% interpolate cube - cube is now stored in Y X Z 
coarseRadDepthCube          = (interp3(Y,X,Z,radDepthCube,Yq,Xq,Zq));
radDepthVcoarse{ctScenNum}  = coarseRadDepthCube(Vcoarse);

end

