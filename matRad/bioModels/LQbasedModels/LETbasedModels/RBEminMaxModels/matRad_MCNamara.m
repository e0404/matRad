classdef matRad_MCNamara < matRad_RBEminMax
    % subclass that implements the MCN model
    % (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4634882/) (accessed on 21/7/2023)
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
        model = 'MCN'
        p0 = 0.99064
        p1 = 0.35605
        p2 = 1.1012
        p3 = -0.0038703

        p0se = 0.014125
        p1se = 0.015038
        p2se = 0.0059972
        p3se = -0.0038703
        requiredQuantities = {'physicalDose', 'LET'}
        possibleRadiationModes = {'protons'}
    end

    methods

        function this = matRad_MCNamara()
            this@matRad_RBEminMax();

        end

        function [RBEmin, RBEmax] = getRBEminMax(this, bixel)

            RBEmax     = this.p0 + ((this.p1 * bixel.LET) ./ bixel.vABratio);
            RBEmin     = this.p2 + (this.p3  * sqrt(bixel.vABratio) .* bixel.LET);

        end

    end
end
