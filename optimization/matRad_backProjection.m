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
    
    % pre-allocation
    d = cell(dij.numOfScenarios,1);
    
    % Calculate dose vector
    if isequal(type,'none')
        
        for i = 1:dij.numOfScenarios
            d{i} = dij.physicalDose{i} * w;
        end
        
    elseif isequal(type,'effect') || isequal(type,'RBExD') 
        
        for i = 1:dij.numOfScenarios
            
            % calculate effect
            linTerm  = dij.mAlphaDose{i} * w;
            quadTerm = dij.mSqrtBetaDose{i} * w;
            e        = linTerm + quadTerm.^2;   

            if ~isequal(type,'RBExD')
                d{i} = e;
            else
                
                % calculate RBX x dose
                scaledEffectSq = (e./dij.bx)+(dij.gamma.^2);
                scaledEffect   = zeros(length(scaledEffectSq),1);
                % compute sqrt(scaledEffect) only for numeric values (not nan) to save time
                [idx,~]           = find(~isnan(scaledEffectSq));
                scaledEffect(idx) = sqrt(scaledEffectSq(idx));
                d{i}              = scaledEffect - dij.gamma;
                
            end
            
       end       
       
    end   
    
    matRad_global_d = d;
    
end

end

