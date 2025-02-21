classdef (Abstract) matRad_LQLETbasedModel < matRad_LQBasedModel
%  This is an Abstract class implementing Linear Quadratic based biological
%  models based on LET calculation. These models require an engine capable
%  of performing LET calculation
%
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

    methods
        function this = matRad_LQLETbasedModel()
            this@matRad_LQBasedModel();
        end

        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernel)
            
            bixel = calcBiologicalQuantitiesForBixel@matRad_LQBasedModel(this,bixel);

            % This LET array is already interpolated on the radiological
            % depths
            
            bixel.LET = kernel.LET;

        end
    end
end