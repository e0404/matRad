classdef matRad_Carabe < matRad_RBEminMax
% subclass that implements the CAR model
% (https://www.tandfonline.com/doi/full/10.1080/09553000601087176?journalCode=irab20)
% (accessed on 21/7/2023)
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

        model = 'CAR';
        p0_CAR = 0.843;
        p1_CAR = 0.154;
        p2_CAR = 2.686;
        p3_CAR = 1.09;
        p4_CAR = 0.006;

        requiredQuantities = {'physicalDose','LET'};
        possibleRadiationModes = {'protons'};

    end

    methods
        function this = matRad_Carabe()
            this@matRad_RBEminMax();
        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)

            RBEmax     = this.p0_CAR + ((this.p1_CAR * this.p2_CAR)./bixel.vABratio ) .* bixel.LET;
            RBEmin     = this.p3_CAR + ((this.p4_CAR * this.p2_CAR)./bixel.vABratio)  .* bixel.LET;

        end
    end
end