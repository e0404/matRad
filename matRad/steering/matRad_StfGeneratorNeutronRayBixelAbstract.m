classdef (Abstract) matRad_StfGeneratorNeutronRayBixelAbstract < matRad_StfGeneratorExternalRayBixelAbstract
% matRad_StfGeneratorNeutronRayBixelAbstract: Abstract Superclass for
% neutron
%   Stf Generators using the ray-bixel mechanism 
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
    
    methods 
        function this = matRad_StfGeneratorNeutronRayBixelAbstract(pln)
            % Constructs ExternalStfGenerator with or without pln
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorExternalRayBixelAbstract(pln);
         end

         function setDefaults(this)
            % Set default values for ExternalStfGenerator
            this.setDefaults@matRad_StfGeneratorExternalRayBixelAbstract();
         end
    end

    methods (Access = protected)
                
        function beam = initBeamData(this,beam)
            beam = this.initBeamData@matRad_StfGeneratorExternalRayBixelAbstract(beam);
            % beam.SCD = this.machine.meta.SCD;
        end

        function beam = setBeamletEnergies(this,beam)
            %Assigns the max photon machine energy to all rays            
            numOfRays = numel(beam.ray);

            for j = numOfRays:-1:1
                beam.ray(j).energy = this.machine.data.energy;
            end
        end

        function beam = initRays(this,beam) 
            %Initializes the geometrical beamlet information for photon bixels (ray corners at isocenter and collimator plane)

            beam = this.initRays@matRad_StfGeneratorExternalRayBixelAbstract(beam);

            rotMat_vectors_T = transpose(matRad_getRotationMatrix(beam.gantryAngle,beam.couchAngle));

            numOfRays = numel(beam.ray);

            %photon ray-target position
            for j = 1:numOfRays
                        beam.ray(j).beamletCornersAtIso = [beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,-beam.bixelWidth/2];...
                            beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,-beam.bixelWidth/2]]*rotMat_vectors_T;
                        % beam.ray(j).rayCorners_SCD = (repmat([0, beam.SCD - beam.SAD, 0],4,1)+ (beam.SCD/beam.SAD) * ...
                        %     [beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                        %     beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,+beam.bixelWidth/2];...
                        %     beam.ray(j).rayPos_bev + [-beam.bixelWidth/2,0,-beam.bixelWidth/2];...
                        %     beam.ray(j).rayPos_bev + [+beam.bixelWidth/2,0,-beam.bixelWidth/2]])*rotMat_vectors_T;
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
            [available,msg] = matRad_StfGeneratorExternalRayBixelAbstract.isAvailable(pln,machine);
           
            if ~available
                return;
            end

            %available = available && isfield(machine.data,'energy') && isscalar(machine.data.energy);
            
            %available = available && isfield(machine.meta,'SCD') && isscalar(machine.meta.SCD);


            if ~available
                msg = 'Your machine file is invalid and does not contain the basic fields required for photon machines!';                
            else
                msg = [];
            end
        end
    end
end

