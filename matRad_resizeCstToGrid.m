function cst = matRad_resizeCstToGrid(cst,vXgridOld,vYgridOld,vZgridOld,vXgridNew,vYgridNew,vZgridNew)
% matRad function to resize the ct to a given resolution
% 
% call
%   cst = matRad_resizeCstToDoseResolution(cst,vXgridOld,vYgridOld,vZgridOld,vXgridNew,vYgridNew,vZgridNew)
%
% input
%   cst:         matRad cst struct
%   vXgridOld:   vector containing old spatial grid points in x [mm] 
%   vYgridOld:   vector containing old spatial grid points in y [mm] 
%   vZgridOld:   vector containing old spatial grid points in z [mm] 
%   vXgridNew:   vector containing new spatial grid points in x [mm] 
%   vYgridNew:   vector containing new spatial grid points in y [mm] 
%   vZgridNew:   vector containing new spatial grid points in z [mm] 
%
% output
%   cst:        updated matRad cst struct containing new linear voxel indices
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:size(cst,1)
   for j = 1:numel(cst{i,4})
      tmpCube              = zeros([numel(vYgridOld) numel(vXgridOld) numel(vZgridOld)]);
      tmpCube(cst{i,4}{j}) = 1;
      cst{i,4}{j}          = find(matRad_interp3(vXgridOld,vYgridOld,vZgridOld, ...
                                                 tmpCube, ...
                                                 vXgridNew,vYgridNew',vZgridNew,'nearest'));
   end
end