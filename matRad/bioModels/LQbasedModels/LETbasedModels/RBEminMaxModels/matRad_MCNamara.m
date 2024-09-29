classdef matRad_MCNamara < matRad_RBEminMax
% subclass that implements the MCN model
% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4634882/) (accessed on 21/7/2023)
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
        model = 'MCN';
        p0_MCN = 0.999064;
        p1_MCN = 0.35605;
        p2_MCN = 1.1012;
        p3_MCN = -0.0038703;
        
        requiredQuantities = {'physicalDose','LET'};
        possibleRadiationModes = {'protons'};
    end

    methods
        
        function this = matRad_MCNamara()
            this@matRad_RBEminMax();
            this.possibleRadiationModes = {'protons'};

        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)

            RBEmax     = this.p0_MCN + ((this.p1_MCN * bixel.LET )./ bixel.vABratio);
            RBEmin     = this.p2_MCN + (this.p3_MCN  * sqrt(bixel.vABratio) .* bixel.LET);

        end

    end
end