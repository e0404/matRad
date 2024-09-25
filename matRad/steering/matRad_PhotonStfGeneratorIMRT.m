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

        function beam = initBeamData(this,beam)
            beam = this.initBeamData@matRad_ExternalStfGeneratorIMRT(beam);
            beam.SCD = this.machine.meta.SCD;
        end

        function beam = initRays(this,beam) 
            %Initializes the geometrical beamlet information for photon bixels (ray corners at isocenter and collimator plane)

            beam = this.initRays@matRad_ExternalStfGeneratorIMRT(beam);

            rotMat_vectors_T = transpose(matRad_getRotationMatrix(beam.gantryAngle,beam.couchAngle));

            numOfRays = numel(beam.ray);

            %photon ray-target position
            for j = 1:numOfRays
                        beam.ray(j).beamletCornersAtIso = [beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,-beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,-beam.bixelWidth/2]]*rotMat_vectors_T;
                        beam.ray(j).rayCorners_SCD = (repmat([0, beam.SCD - beam.SAD, 0],4,1)+ (beam.SCD/beam.SAD) * ...
                            [beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,-beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,-beam.bixelWidth/2]])*rotMat_vectors_T;
            end
        end
        function beam = setBeamletEnergies(this,beam)
            %Assigns the max photon machine energy to all rays
            
            numOfRays = numel(beam.ray);

            for j = numOfRays:-1:1
                beam.ray(j).energy = this.machine.data.energy;
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
