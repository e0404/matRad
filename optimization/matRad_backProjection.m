function [d,d_exp,Omega] = matRad_backProjection(w,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad back projection function to calculate the current dose-,effect- or
% RBExDose- vector based on the dij struct.
% 
% call
%   [d,vOmega] = matRad_backProjection(w,dij,options)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%   options: option struct defining the type of optimization
%
% output
%   d:       dose vector, effect vector or RBExDose vector 
%   d_exp:   expected dose vector
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
global matRad_global_d_exp;
global matRad_global_Omega;

if isequal(w,matRad_global_x)
    
    % get dose from global variable
    d      = matRad_global_d;
    d_exp  = matRad_global_d_exp;
    Omega  = matRad_global_Omega;
    
else
    
    matRad_global_x = w;
    
    % pre-allocation
    d     = cell(options.numOfScen,1);
    d_exp = cell(1,1);
    Omega = cell(size(cst,1),1);
    
    % calculate integral variance vector
    for i = 1:size(cst,1)
       if ~isempty(cst{i,6})                       
             Omega{i}  = ((cst{i,6}(1).penalty/numel(cst{i,4}{1})) * (cst{i,6}(1).mOmega) * w);
       end
    end
    
    
    % Calculate dose vector
    if ~options.bioOpt
       if isequal(options.model,'none')

           d_exp{1} = dij.physicalDoseExp{1} * w;
           
           for i = 1:length(options.ixForOpt)
               d{i} = dij.physicalDose{options.ixForOpt(i)} * w;
           end     
           
       elseif  isequal(options.model,'constRBE')

           d_exp{1} = dij.physicalDoseExp{1} * (w * dij.RBE);
           
           for i = 1:length(options.ixForOpt)
                d{i} =  dij.physicalDose{options.ixForOpt(i)} * (w * dij.RBE);
           end       
       end
    else
        
        % calculated expected effect
        linTermExp  = dij.mAlphaDoseExp{1}    * w;
        quadTermExp = dij.mSqrtBetaDoseExp{1} * w;
        e_exp       = linTermExp + quadTermExp.^2;
        
        if isequal(options.quantityOpt,'effect')
            d_exp{1} = e_exp;
        elseif isequal(options.quantityOpt,'RBExD')
            % calculate expected RBX x dose
            d_exp{1}             = zeros(dij.numOfVoxels,1);                 
            d_exp{1}(dij.ixDose) = sqrt((e_exp(dij.ixDose)./dij.betaX(dij.ixDose))+(dij.gamma(dij.ixDose).^2)) ...
                                  - dij.gamma(dij.ixDose);          
        end
       
        % loop over all scenarios
        for i = 1:length(options.ixForOpt)
            
            % calculate effect
            linTerm  = dij.mAlphaDose{options.ixForOpt(i)}    * w;
            quadTerm = dij.mSqrtBetaDose{options.ixForOpt(i)} * w;
            e        = linTerm + quadTerm.^2;   
            
            if isequal(options.quantityOpt,'effect')
                d{i}     = e;
            elseif isequal(options.quantityOpt,'RBExD')
                % calculate RBX x dose
                d{i}             = zeros(dij.numOfVoxels,1);
                d{i}(dij.ixDose) = sqrt((e(dij.ixDose)./dij.betaX(dij.ixDose))+(dij.gamma(dij.ixDose).^2)) ...
                                    - dij.gamma(dij.ixDose);                 
            else
               error('matRad: Cannot optimze this quantity')
            end
            
        end       
       
    end   
    
    matRad_global_d     = d;
    matRad_global_d_exp = d_exp;
    matRad_global_Omega = Omega;
    
end

end

