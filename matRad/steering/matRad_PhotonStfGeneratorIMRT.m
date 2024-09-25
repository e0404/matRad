classdef matRad_PhotonStfGeneratorIMRT < matRad_ExternalStfGeneratorIMRT

    properties (Constant)
        name = 'Photon IMRT stf Generator';
        shortName = 'photonIMRT';
        possibleRadiationModes = {'photons'};
    end 

    
    
    methods 
        function this = matRad_PhotonStfGeneratorIMRT(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_ExternalStfGeneratorIMRT(pln);

            if isempty(this.radiationMode)
                this.radiationMode = 'photons';
            end
        end            
    end

    methods (Access = protected)        
        function pbMargin = getPbMargin(this)
            pbMargin = this.bixelWidth;
        end
        function rayPos = initializeRayTargetPosition(this,rayPos,rotMat_vectors_T,SAD) 
            %Initializes the geometrical beamlet information for photon bixels (ray corners at isocenter and collimator plane)

            rayPos = this.initializeRayTargetPosition@matRad_ExternalStfGeneratorIMRT(rayPos,rotMat_vectors_T,SAD);

            %photon ray-target position
            for j = 1:rayPos.numOfRays
                        rayPos.ray(j).beamletCornersAtIso = [this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2]]*rotMat_vectors_T;
                        rayPos.ray(j).rayCorners_SCD = (repmat([0, this.machine.meta.SCD - SAD, 0],4,1)+ (this.machine.meta.SCD/SAD) * ...
                            [this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,+rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [-rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2];...
                            this.rayPos(j,:) + [+rayPos.bixelWidth/2,0,-rayPos.bixelWidth/2]])*rotMat_vectors_T;
            end
        end
        function stfElement = setSourceEnergyOnBeam(this,stfElement)
            %Assigns the max photon machine energy to all rays

            for j = stfElement.numOfRays:-1:1
                stfElement.ray(j).energy = this.machine.data.energy;
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
            [available,msg] = matRad_ExternalStfGeneratorIMRT.IsAvailable(pln,machine);

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
                checkModality = any(strcmp(matRad_PhotonStfGeneratorIMRT.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_PhotonStfGeneratorIMRT.possibleRadiationModes, pln.radiationMode));
                
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

            %check required fields
            checkFields = true;

            checkFields = checkFields && isfield(machine.data,'energy') && isscalar(machine.data.energy);
            
            checkFields = checkFields && isfield(machine.meta,'SCD') && isscalar(machine.meta.SCD);


            available = preCheck && checkFields;
        end
    end
end
