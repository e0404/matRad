classdef (Abstract) matRad_RBEminMax < matRad_LQLETbasedModel
%  This is an Abstract class implementing Linear Quadratic based biological
%  models that share the same RBEmin and RBEmax formalism.
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
        function this = matRad_RBEminMax()
            this@matRad_LQLETbasedModel();
        end

        function [alpha,beta] = getAlphaBeta(this,alphaX,betaX,quantityStruct)
            % Get the model specific RBEmin/RBEmax
            [RBEmin, RBEmax] = this.getRBEminMax(quantityStruct);
            
            alpha = RBEmax.*alphaX;
            beta = RBEmin.^2 .* betaX;
        end

        function [bixel] = calcBiologicalQuantitiesForBixel(this,bixel,kernels)
            % This function implement the standard RBEmin/RBEmax formalism.
            % Specific model parameters are computed by the subclass.

            bixel = calcBiologicalQuantitiesForBixel@matRad_LQLETbasedModel(this,bixel,kernels);
 
            % Get the model specific RBEmin/RBEmax
            [RBEmin, RBEmax] = this.getRBEminMax(bixel);

            bixel.alpha = RBEmax.*bixel.vAlphaX;
            bixel.beta = (RBEmin.^2).*bixel.vBetaX;
        end

        function getRBEminMax(~,~,~,~,~)
            matRad_cfg = MatRad_Config.instance();

            matRad_cfg.dispError('Function getRBEminMax needs to be implemented by specific subclass');
        end
    end


end