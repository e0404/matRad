classdef matRad_ParticleStfGeneratorSingleBeamlet < matRad_ParticleStfGeneratorRayBixelAbstract

    properties 
        energy;
    end

    properties (Constant)
        name = 'Particle Single Beamlet';
        shortName = 'particleSingleBixel';
        possibleRadiationModes = {'protons','helium','carbon'};
    end  
    
    methods 
        function this = matRad_ParticleStfGeneratorSingleBeamlet(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_ParticleStfGeneratorRayBixelAbstract(pln);
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

        function beam = setBeamletEnergies(this,beam)
            if isempty(this.energy)
                % Select energy
            else
                [~,ix] = min(abs(this.energy-this.availableEnergies));                
                beam.ray.energy = this.availableEnergies(ix);
            end
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information            
                   
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_ParticleStfGeneratorRayBixelAbstract.IsAvailable(pln,machine);

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
                checkModality = any(strcmp(matRad_ParticleStfGeneratorSingleBeamlet.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_ParticleStfGeneratorSingleBeamlet.possibleRadiationModes, pln.radiationMode));
                
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
