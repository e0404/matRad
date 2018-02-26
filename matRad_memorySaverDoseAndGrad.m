function vParamOut = matRad_memorySaverDoseAndGrad(vParamIn,vLogical,type,dij,jacobVariables)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the ADDITIONAL dose, bixel gradient, or bixel jacobian/struct contributed
% by voxels in the tail if the memorySaverPhoton options is turned on.
% NOTE: the output of this function must be added to the regular dose or
% bixel gradient/jacobian/struct calculated using the matrix multiplication.
%
% call
%   vParamOut = matRad_memorySaverDoseAndGrad(vParamIn,dij,options,jacobVariables)
%
% input
%   vParamIn:    bixel weight OR gradient/jacobian vector
%   dij:         dose influence matrix
%   type:        string, either 'dose' or 'gradient'
%   options:     option struct for optimization
%
% output
%   vParamOut:    dose or gradient/jacobian vector (bixel space)
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


switch type
    
    case {'gradient'}
        
        vParamOut = zeros(dij.totalNumOfBixels,1);
        
    case {'jacobian'}
        
        currConstraints  = jacobVariables.currConstraints;
        jacobLogical     = jacobVariables.jacobLogical;
        numOfConstraints = numel(jacobLogical);
        
        vParamOut = zeros(1,numOfConstraints*nnz(vLogical));
        
    case {'jacobianStruct'}
        
        offset           = jacobVariables.offset;  
        voxInStructure   = jacobVariables.voxInStructure; 
        numOfConstraints = jacobVariables.numOfConstraints;
        vParamOut        = zeros(1,numOfConstraints*nnz(vLogical));
        
    case{'dose'}
        
        vParamOut = zeros(dij.numOfVoxels,1);
        
    otherwise
        
        matRad_dispToConsole('Option must be a string containing either ''gradient'', ''jacobian'', or ''dose''.\n',[],'error');
end

depthOffset = uint32(0);
tailOffset  = uint32(0);
bixelOffset = 1;

for j = 1:dij.totalNumOfRays
    
    if ~vLogical(j)
        continue
    end
    
    depthInd    = depthOffset+(1:uint32(dij.nDepth(j)));
    depthOffset = depthOffset+uint32(dij.nDepth(j));
    
    if strcmp(type,'jacobian')
        indInSparseVec = repmat(bixelOffset,1,numel(currConstraints))...
            +(currConstraints-1)*nnz(vLogical);
    end
    
    for k = depthInd
        
        tailInd    = tailOffset+(1:uint32(dij.nTailPerDepth(k)));
        tailOffset = tailOffset+uint32(dij.nTailPerDepth(k));
        
        voxInd     = dij.ixTail(tailInd);
        
        switch type
            case{'gradient'}
                % vParamIn is g (voxel), vParamOut is g (bixel)
                vParamOut(j) = vParamOut(j) + dij.scaleFactor .* sum(vParamIn(voxInd)) .*dij.bixelDoseTail(k);
            case{'jacobian'}
                % vParamIn is j (voxel), vParamOut is j (bixel)
                vParamOut(indInSparseVec) = vParamOut(indInSparseVec) + dij.scaleFactor .* sum(vParamIn(voxInd,jacobLogical)) .*dij.bixelDoseTail(k);
            case{'jacobianStruct'}
                % vParamIn is j (voxel), vParamOut is j (bixel)
                voxInd = intersect(voxInd,voxInStructure);
                % if there are no voxels in the tail that are part of the structure, move on
                if isempty(voxInd)
                    continue
                end
                vParamOut(offset+bixelOffset) = 1;
            case {'dose'}
                % vParamIn is w, vParamOut is d
                vParamOut(voxInd) = vParamOut(voxInd) + dij.scaleFactor .* vParamIn(j) .*dij.bixelDoseTail(k);
            otherwise
        end
        
    end
    
    bixelOffset = bixelOffset+1;
    
end

