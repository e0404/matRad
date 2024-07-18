classdef (Abstract) matRad_RBEminMax < matRad_LQLETbasedModel
    
    properties

    end

    methods
        function this = matRad_RBEminMax()
            this@matRad_LQLETbasedModel();
            this.requiredQuantities = {'LET'};
        end

        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernels)
            
            bixel = calcBiologicalQuantitiesForBixel@matRad_LQLETbasedModel(this,bixel,kernels);
 
            [RBEmin, RBEmax] = this.getRBEminMax(bixel);

            bixel.alpha = RBEmax.*bixel.vAlphaX;
            bixel.beta = (RBEmin.^2).*bixel.vBetaX;
        end

        function getRBEminMax(~,~,~,~,~)
            matRad_cfg = MatRad_Config.instance();

            matRad_cfg.dispError('Function getRBEminMax needs to be implemented by specific subclass');
        end
    end


end