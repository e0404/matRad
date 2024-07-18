classdef (Abstract) matRad_LQLETbasedModel < matRad_LQBasedModel

    properties

    end

    methods
        function this = matRad_LQLETbasedModel()
            this@matRad_LQBasedModel();
        end

        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernel)
            
            bixel = calcBiologicalQuantitiesForBixel@matRad_LQBasedModel(this,bixel);

            % This LET array is already interpolated on the radiological
            % depths
            
            bixel.LET = kernel.LET;

        end
    end
end