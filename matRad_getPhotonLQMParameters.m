function [ax,bx] = matRad_getPhotonLQMParameters(cst,numVoxel,ctScen,VdoseGrid)
% matRad function to receive the photon LQM reference parameter 
% 
% call
%   cst = matRad_getPhotonLQMParameters(cst,numVoxel,ctScen,VdoseGrid)
%
% input
%   cst:        matRad cst struct
%   numVoxel:   number of voxels of the dose cube
%   ctScen:     CT scenarios for alpha_x and beta_x should be calculated
%   VdoseGrid:  optional linear index vector that allows to specify subindices
%               for which ax and bx will be computed
%
% output
%   ax:         vector containing for each linear voxel index alpha_x
%   bx:         vector containing for each linear voxel index beta_x
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

ax = zeros(numVoxel,ctScen);
bx = zeros(numVoxel,ctScen);

for i = 1:size(cst,1)
   if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
      for ctScen = 1:numel(cst{i,4})
         
         if exist('VdoseGrid','var')
            isInVdoseGrid = ismember(VdoseGrid,cst{i,4}{ctScen});
            ax(VdoseGrid(isInVdoseGrid)) = cst{i,5}.alphaX;
            bx(VdoseGrid(isInVdoseGrid)) = cst{i,5}.betaX;
         else
            ax(cst{i,4}{ctScen},ctScen) = cst{i,5}.alphaX;
            bx(cst{i,4}{ctScen},ctScen) = cst{i,5}.betaX;
         end
         
      end
   end
end




