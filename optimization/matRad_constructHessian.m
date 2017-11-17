function hessianMatrix = matRad_constructHessian(hessianDiag,dij,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to construct the current hessian matrix using the hessian
% diagonal and the dij struct. For linear and quadratic
% objectives/constraints the hessian diagonal is invariant (zero for linear
% functions) and the use of a global variable for the hessian matrix avoids
% the unnecessary and time-consuming matrix-multiplication for these
% functions
% 
% call
%   hessianMatrix = matRad_constructHessian(hessianDiag,dij,options)
%
% input
%   hessianDiag:    hessian diagonal containing second derivative
%   dij:            dose influence matrix
%   options:        option struct defining the type of optimization
%
% output
%   hessianMatrix:  hessian matrix (numOfBixels x numOfBixels) 
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

global matRad_global_hessianDiag;
global matRad_global_hessianMatrix;

if isequal(hessianDiag,matRad_global_hessianDiag)
    
    % get hessian matrix from global variable
    hessianMatrix = matRad_global_hessianMatrix;
    
else
    
    % set global hessian diagonal
    matRad_global_hessianDiag = hessianDiag;
    
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
    matRad_global_hessianMatrix = hessianMatrix;
    
end

end