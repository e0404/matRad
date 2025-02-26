classdef (Abstract) matRad_LQBasedModel < matRad_BiologicalModel
%  This is an Abstract class implementing Linear Quadratic based biological
%  models. These include all biological models exploiting the variable
%  alpha/beta formalism.
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
    
    properties (Constant)
        defaultReportQuantity = 'RBExDose';
    end
    
    properties
        defaultAlphaX = 0.1;        % Tissue defalut alphaX parameter
        defaultBetaX  = 0.05;       % Tissue default betaX parameter
    end

    methods
        function this = matRad_LQBasedModel()
            this@matRad_BiologicalModel;
        end

        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel)

            % This function only initializes the alpha/beta bixel values
            bixel.alpha = NaN*ones(numel(bixel.radDepths),1);
            bixel.beta  = NaN*ones(numel(bixel.radDepths),1);
            
            bixel.vABratio = bixel.vAlphaX ./ bixel.vBetaX;

        end
    end
end