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
end
