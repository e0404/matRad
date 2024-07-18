classdef matRad_KernelBasedLEM < matRad_LQKernelBasedModel

    properties (Constant)
        model = 'LEM';
    end

    methods
        function this = matRad_KernelBasedLEM()
            this@matRad_LQKernelBasedModel();


            this.requiredQuantities = {'alpha', 'beta'};
            this.availableRadiationModalities = {'carbon'};
        end

        
        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernels)
 
            bixel  = calcBiologicalQuantitiesForBixel@matRad_LQBasedModel(this,bixel);
            
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