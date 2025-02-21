classdef matRad_StfGeneratorNeutronSingleBeamlet < matRad_StfGeneratorNeutronRayBixelAbstract
% matRad_StfGeneratorNeutronSingleBeamlet: 
%   Creates a single beamlet for neutrons
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        name = 'Neutron Single Bixel';
        shortName = 'NeutronSingleBixel';
        possibleRadiationModes = {'neutrons'};
    end  
    
    methods 
        function this = matRad_StfGeneratorNeutronSingleBeamlet(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorNeutronRayBixelAbstract(pln);

            if isempty(this.radiationMode)
                this.radiationMode = 'neutrons';
            end
        end            
    end

    methods (Access = protected)  
        function pbMargin = getPbMargin(this)
            pbMargin = 0;
        end
        
        function rayPos = getRayPositionMatrix(this,beam)
            % see superclass for information
            rayPos = [0 0 0];
        end   
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information            
                   
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_StfGeneratorNeutronRayBixelAbstract.isAvailable(pln,machine);

            if ~available
                return;
            else
                available = false;
                msg = [];
            end
    
            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');
    
                %check modality
                checkModality = any(strcmp(matRad_StfGeneratorNeutronSingleBeamlet.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorNeutronSingleBeamlet.possibleRadiationModes, pln.radiationMode));
                
                %Sanity check compatibility
                if checkModality
                    checkModality = strcmp(machine.meta.radiationMode,pln.radiationMode);
                end
    
                preCheck = checkBasic && checkModality;
    
                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            available = preCheck;
        end
    end
end
