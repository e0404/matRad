classdef matRad_HeliumMairani < matRad_RBEminMax

    properties (Constant)

        model = 'HEL';

        p0_HEL = 1.36938e-1; 
        p1_HEL = 9.73154e-3;
        p2_HEL = 1.51998e-2;
    end

    methods
        function this = matRad_HeliumMairani()
            this@matRad_RBEminMax();

            this.availableRadiationModalities = {'helium'};
        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)

          % data-driven RBE parametrization of helium ions
          % https://iopscience.iop.org/article/10.1088/0031-9155/61/2/888
         
          
          % quadratic fit
          %f_Q      = 8.53959e-4 .* bixelLET.^2;
          %RBEmax_Q = 1 + 2.145e-1  + vABratio.^-1 .* f_Q;
          % linear quadratic fit
          %f_LQ      = 2.91783e-1*bixelLET - 9.525e-4*bixelLET.^2;
          %RBEmax_LQ = 1 + ((1.42057e-1 + (vABratio.^-1)) .* f_LQ);
          % linear exponential fit
          %f_LE      = (2.965e-1 * bixelLET) .* exp(-4.90821e-3 * bixelLET);
          %RBEmax_LE = 1 + ((1.5384e-1  + (vABratio.^-1)) .* f_LE);

          f_QE      = (this.p1_HEL * bixel.LET.^2) .* exp(-this.p2_HEL * bixel.LET);
          RBEmax_QE = 1 + ((this.p0_HEL  + (bixel.vABratio.^-1)) .* f_QE);

          % the linear quadratic fit yielded the best fitting result
          RBEmax = RBEmax_QE;
                  
          RBEmin = 1; % no gain in using fitted parameters over a constant value of 1

        end
    end

end