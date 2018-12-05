function d = matRad_backProjection(w,dij,options)
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
    if isequal(options.bioOpt,'none')
        
        for i = 1:options.numOfScenarios
            d{i} = dij.physicalDose{i} * w;
        end
        
    elseif  isequal(options.ID,'protons_const_RBExD')
        
        for i = 1:options.numOfScenarios
             d{i} =  dij.physicalDose{i} * (w * dij.RBE );
        end
        
    elseif (isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LEMIV_RBExD'))
        
        for i = 1:options.numOfScenarios
            
            % calculate effect
            linTerm  = dij.mAlphaDose{i} * w;
            quadTerm = dij.mSqrtBetaDose{i} * w;
            e        = linTerm + quadTerm.^2;   

            if isequal(options.bioOpt,'LEMIV_effect')
                d{i} = e;
            else
                % calculate RBX x dose
                d{i}             = zeros(dij.doseGrid.numOfVoxels,1);
                d{i}(dij.ixDose) = sqrt((e(dij.ixDose)./dij.bx(dij.ixDose))+(dij.gamma(dij.ixDose).^2)) ...
                                    - dij.gamma(dij.ixDose);
               
            end
            
        end       
       
    end   
    
    matRad_global_d = d;
    
end

end

