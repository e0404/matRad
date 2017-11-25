function hessianDiag = matRad_hessianFunc(dij,d_i,prescription,structure,d_ref)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: hessian function for inverse planning supporting 
% squared underdosage, squared overdosage, squared deviation, mean dose 
% objectives, EUD objectives, DVH objectives, (exact) max dose constraint, 
% (exact) min dose constraint, min mean dose constraint, max mean dose 
% constraint, min EUD constraint, max EUD constraint, max DVH constraint,
% min DVH constraint 
%
% call
%   hessianDiag = matRad_hessianFunc(dij,d_i,prescription,structure,d_ref)
%
% input
%   dij:            dose influence matrix
%   d_i:            dose vector
%   prescription:   matRad objective/constraint struct
%   structure:      structure voxels
%   d_ref:          reference dose
%
% output
%   hessianDiag:    hessian diagonal (i.e. second derivative with respect
%                   to dose) of the objectives and constraints. Needs 
%                   subsequent matrix-multiplication (dij-matrix) for 
%                   hessian matrix (see hessianFunWrapper)
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
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

if isequal(prescription.type, 'square underdosing') 

    % underdose : Dose minus prefered dose
    underdose = d_i - d_ref;
    
    % calculate hessian diagonal
    hessianDiag = sparse(structure(underdose<0), ones(nnz(underdose<0),1), (2 * prescription.penalty / numOfVoxels) * ones(nnz(underdose<0),1), dij.numOfVoxels, 1);
    
elseif isequal(prescription.type, 'square overdosing')

    % overdose : Dose minus prefered dose
    overdose = d_i - d_ref;
    
    % calculate hessian diagonal
    hessianDiag = sparse(structure(overdose>0), ones(nnz(overdose>0),1), (2 * prescription.penalty / numOfVoxels) * ones(nnz(overdose>0),1), dij.numOfVoxels, 1);
    
elseif isequal(prescription.type, 'square deviation')
    
    % calculate hessian diagonal
    hessianDiag = sparse(structure, ones(length(structure),1), (2 * prescription.penalty / numOfVoxels) * ones(length(structure),1), dij.numOfVoxels, 1);
    
elseif isequal(prescription.type, 'mean')
    
    % set all-zero hessian diagonal
    hessianDiag = sparse(dij.numOfVoxels,1);
    
elseif isequal(prescription.type, 'EUD')
        
    % get exponent for EUD
    exponent = prescription.EUD;

    % calculate Hessian
    if sum(d_i.^exponent)>0 && exponent ~= 1
        
        % calculate hessian diagonal
        hessianDiag = sparse(structure, ones(length(structure),1), ...
            prescription.penalty*nthroot(1/numOfVoxels,exponent) * ((1-exponent)*sum(d_i.^exponent)^(1/exponent-2)*d_i.^(2*(exponent-1)) + (exponent-1)*sum(d_i.^exponent)^(1/exponent-1)*d_i.^(exponent-2)), ...
            dij.numOfVoxels, 1);    
    else                
        % set all-zero hessian diagonal
        hessianDiag = sparse(dij.numOfVoxels,1);
    end
    
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
    
    % calculate hessian diagonal
    hessianDiag = sparse(structure(voxel_idx), ones(nnz(voxel_idx),1), (2 * prescription.penalty / numOfVoxels) * ones(nnz(voxel_idx),1), dij.numOfVoxels, 1);
    
elseif isequal(prescription.type, 'max EUD constraint') || ...
       isequal(prescription.type, 'min EUD constraint') 

    % exponenent for EUD constraint
    exponent = constraint.EUD;
    
    % calculate hessian diagonal
    hessianDiag = sparse(structure, ones(length(structure),1), ...
        nthroot(1/numOfVoxels,exponent) * ((1-exponent)*sum(d_i.^exponent)^(1/exponent-2)*d_i.^(2*(exponent-1)) + (exponent-1)*sum(d_i.^exponent)^(1/exponent-1)*d_i.^(exponent-2)), ...
        dij.numOfVoxels, 1);    

elseif isequal(prescription.type, 'max dose constraint')
    
    % Approximated max dose constraint not supported
    error('Exact optimization not available for ''%s''!!! Select exact version!', prescription.type)

elseif isequal(prescription.type, 'min dose constraint')
    
    % Approximated min dose constraint not supported
    error('Exact optimization not available for ''%s''!!! Select exact version!', prescription.type)

elseif isequal(prescription.type, 'max mean dose constraint') || ...
       isequal(prescription.type, 'min mean dose constraint') 
    
    % set all-zero hessian diagonal
    hessianDiag = sparse(dij.numOfVoxels,1);
    
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
    
    % exponential term
    exp_term = exp(2*DVHCScaling*(d_i-d_ref));
    
    % calculate hessian diagonal
    hessianDiag = sparse(structure, ones(length(structure),1), (2 * DVHCScaling)^2/numOfVoxels * (exp_term .* (1 - exp_term)) ./ (exp_term + 1).^3, dij.numOfVoxels, 1);    
 
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

elseif isequal(prescription.type, 'max dose constraint (exact)') || ...
       isequal(prescription.type, 'min dose constraint (exact)')
    
    % set all-zero hessian diagonal
    hessianDiag = sparse(dij.numOfVoxels,1);

elseif isequal(prescription.type, 'max dose objective (exact)') || ...
       isequal(prescription.type, 'min dose objective (exact)')
    
    % set all-zero hessian diagonal
    hessianDiag = sparse(dij.numOfVoxels,1);
    
elseif isequal(prescription.type, 'minimax constraint (exact)') || ...
       isequal(prescription.type, 'maximin constraint (exact)')
    
    % set all-zero hessian diagonal
    hessianDiag = sparse(dij.numOfVoxels,1);
end