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

            linTerm  = dij.mAlphaDose{dij.indexforOpt(i)} * w;
            quadTerm = dij.mSqrtBetaDose{dij.indexforOpt(i)} * w;
            e        = linTerm + quadTerm.^2; 
            if isequal(options.quantity,'effect') 
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
       
end   
    
matRad_global_d = d;
    
end


%% reference calculation for the MCNamaraModel
             
%    dp       = dij.physicalDose{dij.indexforOpt(i)} * w;
%    ix       = dp > 0;
%    LETd     = zeros(dij.numOfVoxels,1);
%    LETd(ix) = (dij.mLETDose{dij.indexforOpt(i)}(ix,:)  * w)./dp(ix);
%              
%    sqab     = zeros(dij.numOfVoxels,1);
%    sqab(ix) = sqrt(dij.abX(ix));
%                
%    part1 =  (options.p0 * ((dij.ax .* delta{i})'*dij.physicalDose{1})) + (options.p0 * options.p1 *  ((dij.bx .*delta{i})'*dij.mLETDose{1}));
%    Fac   =  (2*dij.bx .* delta{1}).*((dp * options.p2) -  ( dp * options.p3  .* sqab .* LETd));
%    part2 = (((options.p2*Fac)' * dij.physicalDose{1}) - ((options.p3 * sqab .* Fac)' *dij.mLETDose{1}));            
%    g = g + (part1 + part2)';
                      
