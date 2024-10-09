classdef matRad_HeliumMairani < matRad_RBEminMax
% This class implements the HEL model
% This is a data-driven RBE parametrization of helium ions
% https://iopscience.iop.org/article/10.1088/0031-9155/61/2/888
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

        model = 'HEL';

        p0_HEL = 1.36938e-1; 
        p1_HEL = 9.73154e-3;
        p2_HEL = 1.51998e-2;

        requiredQuantities = {'LET'};
        possibleRadiationModes = {'helium'};
    end

    methods
        function this = matRad_HeliumMairani()
            this@matRad_RBEminMax();
        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)


          f_QE      = (this.p1_HEL * bixel.LET.^2) .* exp(-this.p2_HEL * bixel.LET);
          RBEmax_QE = 1 + ((this.p0_HEL  + (bixel.vABratio.^-1)) .* f_QE);

          % the linear quadratic fit yielded the best fitting result
          RBEmax = RBEmax_QE;
                  
          RBEmin = 1; % no gain in using fitted parameters over a constant value of 1

        end
    end

end