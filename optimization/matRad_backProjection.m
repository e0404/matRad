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
           
             if isequal(options.radMod,'protons')
               
               dp       = dij.physicalDose{dij.indexforOpt(i)} * w;
               ix       = dp > 0;
               LETd     = zeros(dij.numOfVoxels,1);
               RBEmax   = zeros(dij.numOfVoxels,1);
               RBEmin   = zeros(dij.numOfVoxels,1);
               LETd(ix) = (dij.mLETDose{dij.indexforOpt(i)}(ix,:)  * w)./dp(ix);
             
               RBEmax(ix) = options.p0 + ((options.p1 * LETd(ix) )./ dij.abX(ix));           
               RBEmin(ix) = options.p2 - (options.p3  * real(sqrt(dij.abX(ix))) .* LETd(ix)); 
               
             elseif isequal(options.radMod,'carbon')
               % calculate effect
               linTerm  = dij.mAlphaDose{dij.indexforOpt(i)} * w;
               quadTerm = dij.mSqrtBetaDose{dij.indexforOpt(i)} * w;
               e        = linTerm + quadTerm.^2;     
             end
           
            if isequal(options.quantity,'effect') 
               if isequal(options.radMod,'protons')
                  d{i} = dij.ax.*RBEmax.*dp + dij.bx .* RBEmin.^2 .* dp.^2;
               else
                  d{i} = e;
               end
            elseif isequal(options.quantity,'RBExD') 
              % if isequal(options.radMod,'protons')
                  d_ref  = 0.5 .* (sqrt(dij.abX.^2 + (4*dp.*dij.abX.*RBEmax) + (4*dp.^2 .* RBEmin.^2)) - dij.abX);
              % else
                  linTerm  = dij.mAlphaDose{dij.indexforOpt(i)} * w;
                  quadTerm = dij.mSqrtBetaDose{dij.indexforOpt(i)} * w;
                  e        = linTerm + quadTerm.^2;   
                  % calculate RBX x dose
                  scaledEffectSq = (e./dij.bx)+(dij.gamma.^2);
                  scaledEffect   = zeros(length(scaledEffectSq),1);
                  % compute sqrt(scaledEffect) only for numeric values (not nan) to save time
                  [idx,~]           = find(~isnan(scaledEffectSq));
                  scaledEffect(idx) = sqrt(scaledEffectSq(idx));
                  d{i}              = scaledEffect - dij.gamma;
                  
                 if  sum(abs(d_ref -d{1})) > 1e-6
                    error('error');
                 end
              % end
                            
            else
               
               error('not implemented')

            end
            
        end       
       
    end   
    
    matRad_global_d = d;
    
end

end

