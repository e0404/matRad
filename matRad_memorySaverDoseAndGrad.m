function d_or_g_or_j = matRad_memorySaverDoseAndGrad(w_or_g_or_j,dij,option,jacobVariables)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the ADDITIONAL dose, bixel gradient, or bixel jacobian/struct contributed
% by voxels in the tail if the memorySaverPhoton option is turned on.
% NOTE: the output of this function must be added to the regular dose or
% bixel gradient/jacobian/struct calculated using the matrix multiplication.
%
% call
%   output = matRad_memorySaverDoseAndGrad(w_or_g,dij,option)
%
% input
%   w_or_g_or_j:    bixel weight OR gradient/jacobian vector
%   dij:            dose influence matrix
%   options:        string, either 'dose' or 'gradient'
%
% output
%   d_or_g_or_j:    dose vector or gradient/jacobian vector (bixel space)
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

if strcmp(option,'gradient')
    
    d_or_g_or_j = zeros(dij.totalNumOfBixels,1);
elseif strcmp(option,'jacobian')
    
    currConstraints = jacobVariables.currConstraints;
    jacobLogical = jacobVariables.jacobLogical;
    numOfConstraints = numel(jacobLogical);
    
    d_or_g_or_j = zeros(1,numOfConstraints*nnz(dij.optBixel));
elseif strcmp(option,'jacobianStruct')
    
    offset = jacobVariables.offset;
    voxInStructure = jacobVariables.voxInStructure;
    numOfConstraints = jacobVariables.numOfConstraints;
    d_or_g_or_j = zeros(1,numOfConstraints*nnz(dij.optBixel));
elseif strcmp(option,'dose')
    
    d_or_g_or_j = zeros(dij.numOfVoxels,1);
else
    error('Option must be a string containing either ''gradient'', ''jacobian'', or ''dose''.');
end


depthOffset = uint32(0);
tailOffset = uint32(0);
bixelOffset = 1;

for j = 1:dij.totalNumOfRays
    if ~dij.optBixel(j)
        continue
    end
    depthInd = depthOffset+(1:uint32(dij.nDepth(j)));
    depthOffset = depthOffset+uint32(dij.nDepth(j));
    
    if strcmp(option,'jacobian')
        indInSparseVec = repmat(bixelOffset,1,numel(currConstraints))...
            +(currConstraints-1)*nnz(dij.optBixel);
    end
    
    for k = depthInd
        tailInd = tailOffset+(1:uint32(dij.nTailPerDepth(k)));
        tailOffset = tailOffset+uint32(dij.nTailPerDepth(k));
        
        voxInd = dij.ixTail(tailInd);
        
        if strcmp(option,'gradient')
            % w_or_g_or_j is g (voxel), d_or_g_or_j is g (bixel)
            d_or_g_or_j(j) = d_or_g_or_j(j) + dij.scaleFactor .* sum(w_or_g_or_j(voxInd)) .*dij.bixelDoseTail(k);
        elseif strcmp(option,'jacobian')
            % w_or_g_or_j is j (voxel), d_or_g_or_j is j (bixel)
            d_or_g_or_j(indInSparseVec) = d_or_g_or_j(indInSparseVec) + dij.scaleFactor .* sum(w_or_g_or_j(voxInd,jacobLogical)) .*dij.bixelDoseTail(k);
        elseif strcmp(option,'jacobianStruct')
            % w_or_g_or_j is j (voxel), d_or_g_or_j is j (bixel)
            voxInd = intersect(voxInd,voxInStructure);
            % if there are no voxels in the tail that
            % are part of the structure, move on
            if isempty(voxInd)
                continue
            end
            d_or_g_or_j(offset+bixelOffset) = 1;
        elseif strcmp(option,'dose')
            % w_or_g_or_j is w, d_or_g_or_j is d
            d_or_g_or_j(voxInd) = d_or_g_or_j(voxInd) + dij.scaleFactor .* w_or_g_or_j(j) .*dij.bixelDoseTail(k);
        end
    end
    
    bixelOffset = bixelOffset+1;
end




