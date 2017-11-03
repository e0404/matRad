function [hessianMatrix, hessianDiag_new] = matRad_hessianFunc(dij,d_i,prescription,structure,d_ref)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for inverse planning supporting max dose
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint, 
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint 
% 
% call
%   jacobVec = matRad_jacobFunc(d_i,constraint,d_ref)
%
% input
%   d_i:        dose vector
%   constraint: matRad constraint struct
%   d_ref:      reference dose
%
% output
%   jacobVec:  jacobian vector of constraint for differentation with
%              respect to dose. need subsequent differentation for jacobian
%              in beamlet weights (see jacobFunWrapper)
%
% References
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

% number of voxels
numOfVoxels = numel(d_i);

% initialize hessian diagonal
hessianDiag_new = zeros(dij.numOfVoxels,1); %sparse?

if isequal(prescription.type, 'square underdosing') 

    % underdose : Dose minus prefered dose
    underdose = d_i - d_ref;

    % calculate Hessian
    hessianMatrix = sparse(tril(2 * (dij.physicalDose{1}(structure(underdose<0),:)' * dij.physicalDose{1}(structure(underdose<0),:)) * prescription.penalty / numOfVoxels));
    
    % calculate hessian diagonal
    hessianDiag_new(structure(underdose<0)) = 2 * prescription.penalty / numOfVoxels;
    
elseif isequal(prescription.type, 'square overdosing')

    % overdose : Dose minus prefered dose
    overdose = d_i - d_ref;

    % calculate Hessian
    hessianMatrix = sparse(tril(2 * (dij.physicalDose{1}(structure(overdose>0),:)' * dij.physicalDose{1}(structure(overdose>0),:)) * prescription.penalty / numOfVoxels));
    
    % calculate hessian diagonal
    hessianDiag_new(structure(overdose>0)) = 2 * prescription.penalty / numOfVoxels;
    
elseif isequal(prescription.type, 'square deviation')
    
    % calculate Hessian
    hessianMatrix = sparse(tril(2 * (dij.physicalDose{1}(structure,:)' * dij.physicalDose{1}(structure,:)) * prescription.penalty / numOfVoxels));

    % calculate hessian diagonal
    hessianDiag_new(structure) = 2 * prescription.penalty / numOfVoxels;
    
elseif isequal(prescription.type, 'mean')              

    % calculate Hessian
    hessianMatrix = sparse(zeros(dij.totalNumOfBixels));
    
    % calculate hessian diagonal
    % ALL ZERO, NO ADJUSTMENT
    
elseif isequal(prescription.type, 'EUD')
        
    % get exponent for EUD
    exponent = prescription.EUD;

    % calculate Hessian
    if sum(d_i.^exponent)>0 && exponent ~= 1
        
        % calculate Hessian diagonal
        hessianDiag = prescription.penalty*nthroot(1/numOfVoxels,exponent) * ((1-exponent)*sum(d_i.^exponent)^(1/exponent-2)*d_i.^(2*(exponent-1)) + ...
            (exponent-1)*sum(d_i.^exponent)^(1/exponent-1)*d_i.^(exponent-2));
        
        % construct Hessian matrix
        hessianMatrix = sparse(tril(bsxfun(@times, dij.physicalDose{1}(structure,:)', hessianDiag') * dij.physicalDose{1}(structure,:)));
        
        % calculate hessian diagonal
        hessianDiag_new(structure) = prescription.penalty*nthroot(1/numOfVoxels,exponent) * ((1-exponent)*sum(d_i.^exponent)^(1/exponent-2)*d_i.^(2*(exponent-1)) + ...
            (exponent-1)*sum(d_i.^exponent)^(1/exponent-1)*d_i.^(exponent-2));
    else        
        hessianMatrix = sparse(zeros(dij.totalNumOfBixels));
        
        % calculate hessian diagonal
        % ALL ZERO, NO ADJUSTMENT
    end
    
elseif isequal(prescription.type, 'max EUD constraint') || ...
       isequal(prescription.type, 'min EUD constraint') 

    % exponenent for EUD constraint
    exponent = constraint.EUD;

    % calculate Hessian diagonal
    hessianDiag = nthroot(1/numOfVoxels,exponent) * ((1-exponent)*sum(d_i.^exponent)^(1/exponent-2)*d_i.^(2*(exponent-1)) + ...
        (exponent-1)*sum(d_i.^exponent)^(1/exponent-1)*d_i.^(exponent-2));

    % construct Hessian matrix
    hessianMatrix = sparse(tril(bsxfun(@times, dij.physicalDose{1}(structure,:)', hessianDiag') * dij.physicalDose{1}(structure,:)));
    
    % calculate hessian diagonal
    hessianDiag_new(structure) = nthroot(1/numOfVoxels,exponent) * ((1-exponent)*sum(d_i.^exponent)^(1/exponent-2)*d_i.^(2*(exponent-1)) + ...
        (exponent-1)*sum(d_i.^exponent)^(1/exponent-1)*d_i.^(exponent-2));

elseif isequal(prescription.type, 'max dose constraint')
    
    % Approximated max dose constraint not supported
    error('Exact optimization not available for ''%s''!!! Select exact version!', prescription.type)

elseif isequal(prescription.type, 'min dose constraint')
    
    % Approximated min dose constraint not supported
    error('Exact optimization not available for ''%s''!!! Select exact version!', prescription.type)

elseif isequal(prescription.type, 'max mean dose constraint') || ...
       isequal(prescription.type, 'min mean dose constraint') 
   
    % calculate Hessian
    hessianMatrix = sparse(zeros(dij.totalNumOfBixels));
    
    % calculate hessian diagonal
    % ALL ZERO, NO ADJUSTMENT
    
elseif isequal(prescription.type, 'max DVH constraint') || ...
       isequal(prescription.type, 'min DVH constraint')

    d_i_sort = sort(d_i);

    % calculate scaling
    VoxelRatio   = 1;
    NoVoxels     = max(VoxelRatio*numel(d_i),10);
    absDiffsort  = sort(abs(d_ref - d_i_sort));
    deltaDoseMax = absDiffsort(ceil(NoVoxels/2));

    % calclulate DVHC scaling
    ReferenceVal            = 0.01;
    DVHCScaling             = min((log(1/ReferenceVal-1))/(2*deltaDoseMax),250);
    
    % calculate Hessian diagonal
    exp_term = exp(2*DVHCScaling*(d_i-d_ref));    
    hessianDiag = (2 * DVHCScaling)^2/numOfVoxels * (exp_term .* (1 - exp_term)) ./ (exp_term + 1).^3;

    % construct Hessian matrix
    hessianMatrix = sparse(tril(bsxfun(@times, dij.physicalDose{1}(structure,:)', hessianDiag') * dij.physicalDose{1}(structure,:)));
    
    % calculate hessian diagonal
    hessianDiag_new(structure) = (2 * DVHCScaling)^2/numOfVoxels * (exp_term .* (1 - exp_term)) ./ (exp_term + 1).^3;
 
    % alternative constraint calculation 4/4 %               
    % % get reference Volume
    % refVol = cst{j,6}(k).volume/100;
    %  
    % % calc deviation
    % deviation = d_i - d_ref;
    % 
    % % calc d_ref2: V(d_ref2) = refVol
    % d_ref2 = matRad_calcInversDVH(refVol,d_i);
    % 
    % % apply lower and upper dose limits
    % if isequal(cst{j,6}(k).type, 'max DVH constraint')
    %      deviation(d_i < d_ref | d_i > d_ref2) = 0;
    % elseif isequal(cst{j,6}(k).type, 'min DVH constraint')
    %      deviation(d_i > d_ref | d_i < d_ref2) = 0;
    % end
    %   
    % %jacobVec = ones(size(cst{j,4}));             % linear deviation
    % %jacobVec = 2*deviation;                      % square deviation
    % jacobVec = (1/size(cst{j,4},1))*2*deviation; % square deviation with normalization
    % %jacobVec = 4*(deviation).^3;                  % squared square devioation
    % alternative constraint calculation 4/4 %

elseif isequal(prescription.type, 'max DVH objective') ||...
       isequal(prescription.type, 'min DVH objective')

    % get reference Volume
    refVol = prescription.volume/100;

    % calc d_ref2: V(d_ref2) = refVol
    d_ref2 = matRad_calcInversDVH(refVol,d_i);

    % apply lower and upper dose limits
    if isequal(prescription.type, 'max DVH objective')
         voxel_idx = (d_i >= d_ref & d_i <= d_ref2);
    elseif isequal(prescription.type, 'min DVH objective')
         voxel_idx = (d_i <= d_ref & d_i >= d_ref2);
    end
    
    % calculate Hessian
    hessianMatrix = sparse(tril(2 * (dij.physicalDose{1}(structure(voxel_idx),:)' * dij.physicalDose{1}(structure(voxel_idx),:)) * prescription.penalty / numOfVoxels));
    
    % calculate hessian diagonal
    hessianDiag_new(structure(voxel_idx)) = 2 * prescription.penalty / numOfVoxels;
    
elseif isequal(prescription.type, 'max dose constraint (exact)') || ...
       isequal(prescription.type, 'min dose constraint (exact)')
   
    % calculate Hessian
    hessianMatrix = sparse(zeros(dij.totalNumOfBixels));
    
    % calculate hessian diagonal
    % ALL ZERO, NO ADJUSTMENT
end