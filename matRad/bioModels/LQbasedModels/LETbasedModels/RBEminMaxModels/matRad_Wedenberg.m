classdef matRad_Wedenberg < matRad_RBEminMax
% subclass that implements the WED model
% (https://www.ncbi.nlm.nih.gov/pubmed/22909391) (accessed on 21/7/2023)
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
        model  = 'WED';

        p0_WED = 1;
        p1_WED = 0.434;
        p2_WED = 1;

        requiredQuantities = {'physicalDose','LET'};
        possibleRadiationModes = {'protons'};

    end

    methods
        function this = matRad_Wedenberg()
            this@matRad_RBEminMax();
        end

        function [RBEmin, RBEmax] = getRBEminMax(this,bixel)

            RBEmax     = this.p0_WED + ((this.p1_WED * bixel.LET )./ bixel.vABratio);
            RBEmin     = this.p2_WED;
        
        end
           
    end
end