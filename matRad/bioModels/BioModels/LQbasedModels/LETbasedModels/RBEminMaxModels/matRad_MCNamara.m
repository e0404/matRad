classdef matRad_MCNamara < matRad_RBEminMax
% subclass that implements the MCN model
% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4634882/) (accessed on 21/7/2023)    
    
    properties
        
        
    end
    
    properties (Constant)
        model = 'MCN';
        p0_MCN = 0.999064;
        p1_MCN = 0.35605;
        p2_MCN = 1.1012;
        p3_MCN = -0.0038703;

    end

    methods
        
        function this = matRad_MCNamara()
            this@matRad_RBEminMax();

            %this.model = 'MCN';
            this.availableRadiationModalities = {'protons'};

        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)

            RBEmax     = this.p0_MCN + ((this.p1_MCN * bixel.LET )./ bixel.vABratio);
            RBEmin     = this.p2_MCN + (this.p3_MCN  * sqrt(bixel.vABratio) .* bixel.LET);

        end

    end
end