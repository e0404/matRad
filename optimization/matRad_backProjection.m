function d = matRad_backProjection(w,dij,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad back projection function to calculate the current dose-,effect- or
% RBExDose- vector based on the dij struct.
% 
% call
%   d = matRad_backProjection(w,dij,type)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   type: string determing the type of optimization either 'none','effect'
%         or 'RBExD'
%
% output
%   d:    dose vector,effect vector or RBExDose vector 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global matRad_global_x;
global matRad_global_d;

if isequal(w,matRad_global_x)
    
    % get dose from global variable
    d = matRad_global_d;
    
else
    matRad_global_x = w;
    
    % Calculate dose vector
    if isequal(type,'none')
        
        d = dij.physicalDose * w;
        
    elseif isequal(type,'effect') || isequal(type,'RBExD') 
        
        % calculate effect
        linTerm  = dij.mAlphaDose * w;
        quadTerm = dij.mSqrtBetaDose * w;
        e        = linTerm + quadTerm.^2;   

       if ~isequal(type,'RBExD')
           d = e;
       else
           
           % calculate RBX x dose
           ScaledEffect = (e./dij.bx)+(dij.gamma.^2);
           
           % compute sqrt(ScaledEffect) only for numeric values (not nan) to save time
           [idx,~]           = find(~isnan(ScaledEffect));
           ScaledEffect(idx) = sqrt(ScaledEffect(idx));
           d                 = ScaledEffect-dij.gamma;
           
       end       
       
    end   
    matRad_global_d = d;
end

end

