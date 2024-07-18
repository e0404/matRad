classdef matRad_Wedemberg < matRad_RBEminMax
% subclass that implements the WED model
% (https://www.ncbi.nlm.nih.gov/pubmed/22909391) (accessed on 21/7/2023)    
    properties (Constant)
        model  = 'WED';

        p0_WED = 1;
        p1_WED = 0.434;
        p2_WED = 1;

    end

    methods
        function this = matRad_Wedemberg()
            this@matRad_RBEminMax();
            
            this.availableRadiationModalities = {'protons'};
        
        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)

            RBEmax     = this.p0_WED + ((this.p1_WED * bixel.LET )./ bixel.vABratio);
            RBEmin     = this.p2_WED;
        
        end
           
    end
end