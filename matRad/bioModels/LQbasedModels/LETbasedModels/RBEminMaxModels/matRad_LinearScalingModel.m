classdef matRad_LinearScalingModel < matRad_RBEminMax
    % This class implements the Linear Scaling Model
    % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (FITTED for head and neck patients !)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2023-2026 the matRad development team.
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

        model = 'LSM'
        p_lamda_1_1          = 0.008                    % mh:ignore_style (naming scheme intentionally not enforced, kept for backward compatibility)
        p_corrFacEntranceRBE = 0.5    % [kev/mum]        mh:ignore_style
        p_upperLETThreshold  = 30     % [kev/mum]        mh:ignore_style
        p_lowerLETThreshold  = 0.3    % [kev/mum]        mh:ignore_style

        requiredQuantities = {'physicalDose', 'LET'}
        possibleRadiationModes = {'protons', 'helium', 'carbon'}
    end

    methods

        function this = matRad_LinearScalingModel()
            this@matRad_RBEminMax();
        end

        function [RBEmin, RBEmax] = getRBEminMax(this, bixel)

            RBEmax = NaN * ones(size(bixel.vAlphaX));

            ix = this.p_lowerLETThreshold < bixel.LET & bixel.LET < this.p_upperLETThreshold;

            alpha_0 = bixel.vAlphaX - (this.p_lamda_1_1 * this.p_corrFacEntranceRBE);

            RBEmax(ix)  = alpha_0(ix) + this.p_lamda_1_1 * bixel.LET(ix);

            if sum(ix) < length(bixel.LET)
                ixUpper = bixel.LET >= this.p_upperLETThreshold;
                ixLower = bixel.LET <= this.p_lowerLETThreshold;

                RBEmax(ixUpper) = alpha_0(ixUpper) + this.p_lamda_1_1 * this.p_upperLETThreshold;
                RBEmax(ixLower) = alpha_0(ixLower) + this.p_lamda_1_1 * this.p_lowerLETThreshold;
            end

            RBEmax = RBEmax ./ bixel.vAlphaX;
            RBEmin = 1;

        end

    end
end
