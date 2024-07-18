classdef (Abstract) matRad_BiologicalModel < handle
%  matRad_BiologicalModel
%  This is an abstract interface class to define Biological Models for use in
%  dose calculation and plan optimization.
%  Subclasses should at least implement the methods:
% 
%  calcBiologicalQuantitiesForBixel()      to implement the specific model calculation algorithm
%
%  Additional methods that can be implemented are:
% 
%   getTissueInformation()                  to collect meta information about the tissue defintion and paramters                                        
%   assignBioModelPropertiesFromEngine()    to translate user defined paramters directly to the biological model subclass
% 
% 
% All subclasses should also declare the  properties:
%
%   'availableRadiationModalities'         to specify the radiation modalities to which the model validity is limited
%   'requiredQuantities'                   to check the availability of information stored in the provided machine file
%
% constructor (Abstract)
%   matRad_BiologicalModel()
%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        requiredQuantities;                % kernels in base data needed for the alpha/beta calculation
        availableRadiationModalities;      % radiation modalitites compatible with the model
    end

    properties (Abstract, Constant)
        model;
    end

    methods
        function this = matRad_BiologicaModel()
            
        end

        function bixel = calcBiologicalQuantitiesForBixel(this)
            
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Function: calcBiologicalQuantitiesForBixel should be implemented by the model subclass!');
        end

        function assignBioModelPropertiesFromEngine(this, pln)
            
            % This function can be implemented by the specific subclasses
            % to assign model-specific user defined paramters

        end
        
    end

    methods %(Static)
        
        function [vTissueIndex] = getTissueInformation(this,~,~,~,vAlphaX,~,~,~) %(machine,cst,dij,vAlphaX,vBetaX,VdoseGrid, VdoseGridScenIdx)
            
            % This is the default, should be masked by the specific model
            % subclass if needed

            for s=1:numel(vAlphaX)
                vTissueIndex{s} = zeros(size(vAlphaX{s}));
            end


        
        end


    end

end