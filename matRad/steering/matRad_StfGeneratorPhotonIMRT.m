classdef matRad_StfGeneratorPhotonIMRT < matRad_StfGeneratorPhotonRayBixelAbstract

    properties (Constant)
        name = 'Photon IMRT stf Generator';
        shortName = 'PhotonIMRT';
        possibleRadiationModes = {'photons'};
    end 

    
    
    methods 
        function this = matRad_StfGeneratorPhotonIMRT(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorPhotonRayBixelAbstract(pln);

            if isempty(this.radiationMode)
                this.radiationMode = 'photons';
            end
        end            
    end

    methods (Access = protected)        
        function pbMargin = getPbMargin(this)
            pbMargin = this.bixelWidth;
        end        
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information            
                   
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_StfGeneratorPhotonRayBixelAbstract.isAvailable(pln,machine);

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
                checkModality = any(strcmp(matRad_StfGeneratorPhotonIMRT.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorPhotonIMRT.possibleRadiationModes, pln.radiationMode));
                
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
