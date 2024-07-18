classdef matRad_Carabe < matRad_RBEminMax
% subclass that implements the CAR model
% (https://www.tandfonline.com/doi/full/10.1080/09553000601087176?journalCode=irab20)
% (accessed on 21/7/2023)
    properties (Constant)

        model = 'CAR';
        p0_CAR = 0.843;
        p1_CAR = 0.154;
        p2_CAR = 2.686;
        p3_CAR = 1.09;
        p4_CAR = 0.006;

    end

    methods
        function this = matRad_Carabe()
            this@matRad_RBEminMax();
            this.availableRadiationModalities = {'protons'};
        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)

            RBEmax     = this.p0_CAR + ((this.p1_CAR * this.p2_CAR)./bixel.vABratio ) .* bixel.LET;
            RBEmin     = this.p3_CAR + ((this.p4_CAR * this.p2_CAR)./bixel.vABratio)  .* bixel.LET;

        end
    end
end