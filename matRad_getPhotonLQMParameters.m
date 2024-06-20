function [ax,bx] = matRad_getPhotonLQMParameters(cst,numVoxel,VdoseGrid)
% matRad function to receive the photon LQM reference parameter
%
% call
%   [ax,bx] = matRad_getPhotonLQMParameters(cst,numVoxel,ctScen,VdoseGrid)
%
% input
%   cst:            matRad cst struct
%   numVoxel:       number of voxels of the dose cube
%   VdoseGrid:      optional linear index vector that allows to specify subindices
%                   for which ax and bx will be computed
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

numOfCtScen = unique(cellfun(@numel,cst(:,4)));

if numel(numOfCtScen) > 1
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Inconsinstent number of ct scnearios in cst!');
end

ax = cell(numOfCtScen,1);
bx = cell(numOfCtScen,1);

[ax{:}] = deal(zeros(numVoxel,1));
[bx{:}] = deal(zeros(numVoxel,1));

for i = 1:size(cst,1)
    if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
        for s = 1:numOfCtScen

            if exist('VdoseGrid','var')
                if iscell(VdoseGrid)
                    isInVdoseGrid = ismember(VdoseGrid{s},cst{i,4}{s});
                else
                    isInVdoseGrid = ismember(VdoseGrid,cst{i,4}{s});
                end
                ax{s}(VdoseGrid(isInVdoseGrid)) = cst{i,5}.alphaX;
                bx{s}(VdoseGrid(isInVdoseGrid)) = cst{i,5}.betaX;
            else
                ax{s}(cst{i,4}{s}) = cst{i,5}.alphaX;
                bx{s}(cst{i,4}{s}) = cst{i,5}.betaX;
            end

        end
    end
end