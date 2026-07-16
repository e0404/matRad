classdef matRad_MKM < matRad_LQZStarBasedModel
% subclass that implements the MKM model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    properties (Constant)
        model = 'MKM';        
        requiredQuantities = {'physicalDose'};
        possibleRadiationModes = {'protons','helium','carbon','oxygen'};
        kernelQuantities = {'zs'};
    end

    methods
    


      function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernels)
                
          bixel  = calcBiologicalQuantitiesForBixel@matRad_LQZStarBasedModel(this,bixel,kernels);                            
          bixel.alpha = bixel.vAlphaX + bixel.vBetaX.*kernels.zs;
          bixel.beta  = bixel.vBetaX;
      end


    end
end
