classdef matRad_StfGeneratorParticleVHEE < matRad_StfGeneratorParticleRayBixelAbstract
% matRad_StfGeneratorParticleVHEE: VHEE Steering Geometry Setup (stf)
%   Creates the stf data structure containing the steering information /
%   field geometry for VHEE plans on a regular lateral spot grid by using a
%   single, manually defined energy setting.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2025 the matRad development team.
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
        name = 'VHEE Geometry Generator';
        shortName = 'VHEE';
        possibleRadiationModes = {'VHEE'};
    end 

    properties
        energy;
    end

    methods 
        function this = matRad_StfGeneratorParticleVHEE(pln)
            if nargin < 1
                pln = [];
            end

            this@matRad_StfGeneratorParticleRayBixelAbstract(pln);
         end
    end


    methods (Access = protected)

        function beam = initBeamData(this, beam)
            % Initialize beam data from the superclass, then handle single-energy logic
            beam = initBeamData@matRad_StfGeneratorParticleRayBixelAbstract(this, beam);
            % If user doesn't specify an energy, default to 200 MeV
            if isempty(this.energy)
                beam.VHEEenergy = 200;  % Default
            else
                beam.VHEEenergy = this.energy;
            end
            this.energy = beam.VHEEenergy;
            % Optional: check if that energy is in the machine data
            if isfield(this.machine.data,'energies') && ~isempty(this.machine.data.energies)
                if ~ismember(beam.VHEEenergy, this.machine.data.energies)
                    error(['The specified VHEE energy (',num2str(beam.VHEEenergy), ...
                           ' MeV) is not found in machine.data.energies!']);
                end
            end
        end

        function beam = setBeamletEnergies(~,beam) 
            %Assigns defined particle machine energy from plan to all rays
           
            beam.numOfBixelsPerRay = zeros(1,beam.numOfRays);
            beam.numOfRays = numel(beam.ray);
            for j = beam.numOfRays:-1:1

                %fix the energy for all rays for VHEE
                beam.ray(j).focusIx = 1;
                
                beam.ray(j).energy = beam.VHEEenergy;
                beam.ray(j).rangeShifter.ID = 0;
                beam.ray(j).rangeShifter.eqThickness = 0;
                beam.ray(j).rangeShifter.sourceRashiDistance = 0;
              
                beam.numOfBixelsPerRay(j) = 1;
            end
        end

        function  beam = finalizeBeam(this,beam)
            for j = beam.numOfRays:-1:1
                if isempty(beam.ray(j).energy)
                    beam.ray(j) = [];
                    beam.numOfBixelsPerRay(j) = [];
                    beam.numOfRays = beam.numOfRays - 1;
                end
            end
            beam = this.finalizeBeam@matRad_StfGeneratorParticleRayBixelAbstract(beam);
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information            
                   
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_StfGeneratorParticleRayBixelAbstract.isAvailable(pln,machine);

            if ~available
                return;
            else
                available = false;
                msg = [];
            end
    
            %checkBasic
            try    
                %check modality
                checkModality = any(strcmp(matRad_StfGeneratorParticleVHEE.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorParticleVHEE.possibleRadiationModes, pln.radiationMode));
                
                %Sanity check compatibility
                if checkModality
                    checkModality = strcmp(machine.meta.radiationMode,pln.radiationMode);
                end
    
                preCheck = checkModality;
    
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
