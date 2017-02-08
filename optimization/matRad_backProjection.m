function d = matRad_backProjection(w,dij,options)
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

global matRad_global_x;
global matRad_global_d;

if isequal(w,matRad_global_x)
    
    % get dose from global variable
    d = matRad_global_d;
    
else
    
    matRad_global_x = w;
    
    % pre-allocation
    d = cell(options.numOfScenarios,1);
    
    % Calculate dose vector
    if isequal(options.quantity,'physicalDose')
        for i = 1: length(dij.indexforOpt)
            d{i} = dij.physicalDose{dij.indexforOpt(i)} * w;
        end
        
    elseif  isequal(options.type,'const_RBExD')
        
        for i = 1:length(dij.indexforOpt)
             d{i} =  dij.physicalDose{dij.indexforOpt(i)} * (w * dij.RBE );
        end
        
    elseif options.bioOpt
        
        for i = 1:length(dij.indexforOpt)
            
            % calculate effect
            linTerm  = dij.mAlphaDose{dij.indexforOpt(i)} * w;
            quadTerm = dij.mSqrtBetaDose{dij.indexforOpt(i)} * w;
            e        = linTerm + quadTerm.^2;   
            
            if isequal(options.quantity,'effect') 
                d{i} = e;
            elseif isequal(options.quantity,'RBExD') 
                % calculate RBX x dose
                scaledEffectSq = (e./dij.bx)+(dij.gamma.^2);
                scaledEffect   = zeros(length(scaledEffectSq),1);
                % compute sqrt(scaledEffect) only for numeric values (not nan) to save time
                [idx,~]           = find(~isnan(scaledEffectSq));
                scaledEffect(idx) = sqrt(scaledEffectSq(idx));
                d{i}              = scaledEffect - dij.gamma;
                
%                 ab = dij.ax./dij.bx;
%                 dp = dij.physicalDose{dij.indexforOpt(i)} * w;
%                 let = (dij.mLETDose{dij.indexforOpt(i)} * w)./dp;
%                 RBEmax = 0.999064 + ((0.35605 ./ab).*let);
%                 RBEmin = 1.1012-0.0038703 * sqrt(ab).*let;
%                 
%                 d{i}  = 1./(2) .* (sqrt(ab.^2 + (4*dp.*ab.*RBEmax) + (4*dp.^2 .* RBEmin.^2)) - ab);
%                 
%                 vDiff = d2-d{i};
            else
               
               error('not implemented')

            end
            
        end       
       
    end   
    
    matRad_global_d = d;
    
end

end

