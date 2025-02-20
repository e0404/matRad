classdef matRad_StfGeneratorVHEE < matRad_StfGeneratorParticleRayBixelAbstract
% matRad_ParticleStfGenerator: Abstract Superclass for Steering information 
%   generators. Steering information is used to guide the dose calculation
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
        name = 'VHEE stf Generator';
        shortName = 'ParticleVHEE';
        possibleRadiationModes = {'VHEE'};
    end 

    properties
        longitudinalSpotSpacing;
    end

    methods 
        function this = matRad_StfGeneratorVHEE(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorParticleRayBixelAbstract(pln);
         end
    end


    methods (Access = protected)

        function beam = initBeamData(this,beam)
            beam = this.initBeamData@matRad_StfGeneratorParticleRayBixelAbstract(beam);
            beam.longitudinalSpotSpacing = this.longitudinalSpotSpacing;
        end
        
        function beam = setBeamletEnergies(this,beam) 
            %Assigns the max particle machine energy layers to all rays

            isoCenterInCubeCoords = matRad_world2cubeCoords(beam.isoCenter,this.ct);

            if isfield(this.machine.meta,'LUT_bxWidthminFWHM')
                LUTspotSize = this.machine.meta.LUT_bxWidthminFWHM;
            else
                LUTspotSize = this.machine.meta.LUTspotSize;
            end
            

            beam.numOfBixelsPerRay = zeros(1,beam.numOfRays);

            for j = beam.numOfRays:-1:1

                for shiftScen = 1:this.multScen.totNumShiftScen
                        % ray tracing necessary to determine depth of the target
                        [alphas,l{shiftScen},rho{shiftScen},d12,~] = matRad_siddonRayTracer(isoCenterInCubeCoords + this.multScen.isoShift(shiftScen,:), ...
                            this.ct.resolution, ...
                            beam.sourcePoint, ...
                            beam.ray(j).targetPoint, ...
                            [this.ct.cube {this.voiTarget}]);

                        %Used for generic range-shifter placement
                        ctEntryPoint = alphas(1) * d12;
                end

                % target hit
                rhoVOITarget = [];
                for shiftScen = 1:this.multScen.totNumShiftScen
                    rhoVOITarget = [rhoVOITarget, rho{shiftScen}{end}];
                end

                if any(rhoVOITarget)
                    counter = 0;

                    %Here we iterate through scenarios to check the required
                    %energies w.r.t lateral position.
                    %TODO: iterate over the linear scenario mask instead?

     

                    beam.ray(j).energy = this.pln.VHEE_energy;
                    beam.ray(j).rangeShifter = [];

                    % Save energies in stf struct
                    for k = 1:numel(targetEntry)

                        %If we need lower energies than available, consider
                        %range shifter (if asked for)
                        if any(targetEntry < min(this.availablePeakPos)) && this.useRangeShifter
                            %Get Energies to use with range shifter to fill up
                            %non-reachable low-range spots
                            raShiEnergies = this.availableEnergies(this.availablePeakPosRaShi >= targetEntry(k) & min(this.availablePeakPos) > this.availablePeakPosRaShi);

                            raShi.ID = 1;
                            raShi.eqThickness = rangeShifterEqD;
                            raShi.sourceRashiDistance = round(ctEntryPoint - 2*rangeShifterEqD,-1); %place a little away from entry, round to cms to reduce number of unique settings

                            beam.ray(j).energy = [beam.ray(j).energy raShiEnergies];
                            beam.ray(j).rangeShifter = [beam.ray(j).rangeShifter repmat(raShi,1,length(raShiEnergies))];
                        end

                        %Normal placement without rangeshifter
                        newEnergies = this.availableEnergies(this.availablePeakPos>=targetEntry(k)&this.availablePeakPos<=targetExit(k));


                        beam.ray(j).energy = [beam.ray(j).energy newEnergies];


                        raShi.ID = 0;
                        raShi.eqThickness = 0;
                        raShi.sourceRashiDistance = 0;
                        beam.ray(j).rangeShifter = [beam.ray(j).rangeShifter repmat(raShi,1,length(newEnergies))];
                    end


                    targetEntry = [];
                    targetExit = [];


                    % book keeping & calculate focus index
                    beam.numOfBixelsPerRay(j) = numel([beam.ray(j).energy]);
                    currentMinimumFWHM = matRad_interp1(LUTspotSize(1,:)',...
                        LUTspotSize(2,:)',...
                        this.bixelWidth, ...
                        LUTspotSize(2,end));
                    focusIx  =  ones(beam.numOfBixelsPerRay(j),1);
                    [~, vEnergyIx] = min(abs(bsxfun(@minus,[this.machine.data.energy]',...
                        repmat(beam.ray(j).energy,length([this.machine.data]),1))));

                    % get for each spot the focus index
                    for k = 1:beam.numOfBixelsPerRay(j)
                        focusIx(k) = find(this.machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                    end

                    beam.ray(j).focusIx = focusIx';

                    %Get machine bounds
                    numParticlesPerMU = 1e6*ones(1,beam.numOfBixelsPerRay(j));
                    minMU = zeros(1,beam.numOfBixelsPerRay(j));
                    maxMU = Inf(1,beam.numOfBixelsPerRay(j));
                    for k = 1:beam.numOfBixelsPerRay(j)
                        if isfield(this.machine.data(vEnergyIx(k)),'MUdata')
                            MUdata = this.machine.data(vEnergyIx(k)).MUdata;
                            if isfield(MUdata,'numParticlesPerMU')
                                numParticlesPerMU(k) = MUdata.numParticlesPerMU;
                            end

                            if isfield(MUdata,'minMU')
                                minMU(k) = MUdata.minMU;
                            end

                            if isfield(MUdata,'maxMU')
                                maxMU(k) = MUdata.maxMU;
                            end
                        end
                    end


                else % target not hit
                    beam.ray(j)               = [];
                    beam.numOfBixelsPerRay(j) = [];
                end
            end
            beam.numOfRays = numel(beam.ray);
        end

        function  beam = finalizeBeam(this,beam)

            for j = beam.numOfRays:-1:1
                for k = beam.numOfBixelsPerRay(j):-1:1
                    maskEnergy = beam.ray(j).energy(k) == this.availableEnergies;
                    if ~useEnergyBool(maskEnergy)
                        beam.ray(j).energy(k)         = [];
                        beam.ray(j).focusIx(k)        = [];
                        beam.ray(j).rangeShifter(k)   = [];
                        beam.numOfBixelsPerRay(j) = beam.numOfBixelsPerRay(j) - 1;
                    end
                end
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
                checkModality = any(strcmp(matRad_StfGeneratorVHEE.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorVHEE.possibleRadiationModes, pln.radiationMode));
                
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

