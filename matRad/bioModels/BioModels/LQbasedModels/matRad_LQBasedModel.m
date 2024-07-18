classdef (Abstract) matRad_LQBasedModel < matRad_BiologicalModel

    properties
        defaultAlphaX = 0.1;
        defaultBetaX  = 0.05;
    end

    methods
        function this = matRad_LQBasedModel()
            this@matRad_BiologicalModel;
        end

        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel)

            % This function only initializes the alpha/beta bixel values
            bixel.alpha = NaN*ones(numel(bixel.radDepths),1);
            bixel.beta  = NaN*ones(numel(bixel.radDepths),1);
            
            bixel.vABratio = bixel.vAlphaX ./ bixel.vBetaX;

        end
    end
end