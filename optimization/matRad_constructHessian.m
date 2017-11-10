function hessianMatrix = matRad_constructHessian(hessianDiag,dij,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad back projection function to calculate the current dose-,effect- or
% RBExDose- vector based on the dij struct.
% 
% call
%   d = matRad_backProjection(w,dij,options)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   options: option struct defining the type of optimization
%
% output
%   d:       dose vector, effect vector or RBExDose vector 
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

global matRad_global_hessian_diag;
global matRad_global_hessian_matrix;

if isequal(hessianDiag,matRad_global_hessian_diag)
    
    % get hessian from global variable
    hessianMatrix = matRad_global_hessian_matrix;
    
else
    
    % set global hessian diagonal
    matRad_global_hessian_diag = hessianDiag;
    
    % construct hessian matrix
    if isequal(options.bioOpt,'none')
        
        voxel_idx = any(hessianDiag,2) & any(dij.physicalDose{1},2);
        hessianMatrix = sparse(tril(bsxfun(@times, dij.physicalDose{1}(voxel_idx,:)', hessianDiag(voxel_idx)') * dij.physicalDose{1}(voxel_idx,:)));
        
    elseif  isequal(options.ID,'protons_const_RBExD')
        
        error('biological optimization not supported yet for exact optimization!!!')
        
    elseif (isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LEMIV_RBExD'))
        
        error('biological optimization not supported yet for exact optimization!!!')  
       
    end   
    
    % set global hessian matrix
    matRad_global_hessian_matrix = hessianMatrix;
    
end

end