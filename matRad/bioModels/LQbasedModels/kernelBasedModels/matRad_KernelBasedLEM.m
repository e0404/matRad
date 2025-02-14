classdef matRad_KernelBasedLEM < matRad_LQKernelBasedModel
% This class specifically implements the kernel-based LEMIV model
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
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
        model = 'LEM';
        
        requiredQuantities = {'alpha','beta','physicalDose'};
        kernelQuantities = {'alpha','beta'};

        possibleRadiationModes = {'protons','helium','carbon'};
    end

    methods
        function this = matRad_KernelBasedLEM()
            this@matRad_LQKernelBasedModel();
        end

        
        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernels)
            % This function assignis the interpolated kernels to the
            % correct tissue class
            
            bixel  = calcBiologicalQuantitiesForBixel@matRad_LQKernelBasedModel(this,bixel);
            
            numOfTissueClass = size(bixel.baseData.alpha,2);
            
            for i = 1:numOfTissueClass
                mask = bixel.vTissueIndex == i;
                if any(mask)
                    bixel.alpha(mask) = kernels.alpha(mask,i);
                    bixel.beta(mask)  = kernels.beta(mask,i);
                end
            end       
        end
    end
end