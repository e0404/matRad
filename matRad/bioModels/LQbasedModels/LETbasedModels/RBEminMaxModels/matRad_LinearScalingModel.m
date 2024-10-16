classdef matRad_LinearScalingModel < matRad_RBEminMax
% This class implements the Linear Scaling Model
% according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (FITTED for head and neck patients !)
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

        model = 'LSM';
        p_lamda_1_1          = 0.008;
        p_corrFacEntranceRBE = 0.5;   %[kev/mum]
        p_upperLETThreshold  = 30;    %[kev/mum]
        p_lowerLETThreshold  = 0.3;   %[kev/mum]

        requiredQuantities = {'physicalDose','LET'};
        possibleRadiationModes = {'protons','helium','carbon'};
    end

    methods
        function this = matRad_LinearScalingModel()
            this@matRad_RBEminMax();
        end

        function [RBEmin,RBEmax] = getRBEminMax(this,bixel)

            RBEmax = NaN*ones(size(bixel.vAlphaX));
            
            ix = this.p_lowerLETThreshold < bixel.LET < this.p_upperLETThreshold;
                  
            alpha_0 = bixel.vAlphaX - (this.p_lamda_1_1 * this.p_corrFacEntranceRBE);
            
            RBEmax(ix)  = alpha_0(ix) + this.p_lamda_1_1 * bixel.LET;
            
             if sum(ix) < length(bixel.LET)
                 RBEmax(bixel.LET > this.p_upperLETThreshold) =  alpha_0(bixel.LET > this.p_upperLETThreshold) + this.p_lamda_1_1 * this.p_upperLETThreshold;
                 RBEmax(bixel.LET < this.p_lowerLETThreshold) =  alpha_0(bixel.LET < this.p_lowerLETThreshold) + this.p_lamda_1_1 * this.p_lowerLETThreshold;
             end

             RBEmax = RBEmax./bixel.vAlphaX(ix);
             RBEmin = 1;
            
        end

    end
end