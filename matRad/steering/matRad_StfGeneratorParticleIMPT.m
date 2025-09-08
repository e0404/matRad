classdef matRad_StfGeneratorParticleIMPT < matRad_StfGeneratorParticleRayBixelAbstract
% matRad_StfGeneratorParticleIMPT: IMPT Steering Geometry Setup (stf)
%   Creates the stf data structure containing the steering information /
%   field geometry for standard IMPT plans on a regular lateral spot grid
%   with modulation in depth through energy layers
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
        name = 'Particle IMPT stf Generator';
        shortName = 'ParticleIMPT';
        possibleRadiationModes = {'protons','helium','carbon'};
        airOffsetCorrection = true;
    end 

    properties
        longitudinalSpotSpacing;
    end

    methods 
        function this = matRad_StfGeneratorParticleIMPT(pln)
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
            matRad_cfg = MatRad_Config.instance();

            isoCenterInCubeCoords = matRad_world2cubeCoords(beam.isoCenter,this.ct);

            if isfield(this.machine.meta,'LUT_bxWidthminFWHM')
                LUTspotSize = this.machine.meta.LUT_bxWidthminFWHM;
            else
                LUTspotSize = this.machine.meta.LUTspotSize;
            end

            %Air Offset Correction
            if this.airOffsetCorrection
                if ~isfield(this.machine.meta, 'fitAirOffset')
                    this.machine.meta.fitAirOffset = 0; %By default we assume that the base data was fitted to a phantom with surface at isocenter
                    matRad_cfg.dispDebug('Asked for correction of Base Data Air Offset, but no value found. Using default value of %f mm.\n',this.machine.meta.fitAirOffset);
                end
            else
                this.machine.meta.fitAirOffset = 0;
            end

            beam.numOfBixelsPerRay = zeros(1,beam.numOfRays);

            for j = beam.numOfRays:-1:1
                
                ctEntryPoint = zeros(this.multScen.totNumShiftScen,1);

                radDepthOffset = zeros(this.multScen.totNumShiftScen,1);

                for shiftScen = 1:this.multScen.totNumShiftScen
                        % ray tracing necessary to determine depth of the target
                        [alphas,l{shiftScen},rho{shiftScen},d12,~] = matRad_siddonRayTracer(isoCenterInCubeCoords + this.multScen.isoShift(shiftScen,:), ...
                            this.ct.resolution, ...
                            beam.sourcePoint, ...
                            beam.ray(j).targetPoint, ...
                            [this.ct.cube {this.voiTarget}]);

                        %Used for generic range-shifter placement
                        ctEntryPoint(shiftScen) = alphas(1) * d12;

                        if this.airOffsetCorrection
                            nozzleToSkin = (ctEntryPoint(shiftScen) + this.machine.meta.BAMStoIsoDist) - this.machine.meta.SAD;
                            radDepthOffset(shiftScen) = 0.0011 * (nozzleToSkin - this.machine.meta.fitAirOffset);
                        else
                            radDepthOffset(shiftScen) = 0;
                        end
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
                    for ctScen = 1:this.multScen.numOfCtScen
                        for shiftScen = 1:this.multScen.totNumShiftScen
                            for rangeShiftScen = 1:this.multScen.totNumRangeScen
                                if this.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                                    counter = counter+1;

                                    % compute radiological depths
                                    % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                                    rSP = l{shiftScen} .* rho{shiftScen}{ctScen} ;
                                    radDepths = cumsum(rSP) - 0.5*rSP + radDepthOffset(shiftScen);

                                    if this.multScen.relRangeShift(rangeShiftScen) ~= 0 || this.multScen.absRangeShift(rangeShiftScen) ~= 0
                                        radDepths = radDepths +...                                                      % original cube
                                            rho{shiftScen}{ctScen}*this.multScen.relRangeShift(rangeShiftScen) +...     % rel range shift
                                            this.multScen.absRangeShift(rangeShiftScen);                                % absolute range shift
                                        radDepths(radDepths < 0) = 0;
                                    end

                                    % find target entry & exit
                                    diff_voi    = [diff([rho{shiftScen}{end}])];
                                    entryIx = find(diff_voi == 1);
                                    exitIx = find(diff_voi == -1);

                                    %We approximate the interface using the rad depth between the last voxel before and the first voxel after the interface 
                                    % This captures the case that the first relevant voxel is a target voxel
                                    targetEntry(counter,1:length(entryIx))  = (radDepths(entryIx) + radDepths(entryIx+1)) ./ 2;
                                    targetExit(counter,1:length(exitIx))    = (radDepths(exitIx) + radDepths(exitIx+1))   ./ 2;
                                end
                            end
                        end
                    end

                    targetEntry(targetEntry == 0) = NaN;
                    targetExit(targetExit == 0)   = NaN;

                    targetEntry = min(targetEntry);
                    targetExit  = max(targetExit);

                    %check that each energy appears only once in stf
                    if(numel(targetEntry)>1)
                        m = numel(targetEntry);
                        while(m>1)
                            if(targetEntry(m) < targetExit(m-1))
                                targetExit(m-1) = max(targetExit(m-1:m));
                                targetExit(m)=[];
                                targetEntry(m-1) = min(targetEntry(m-1:m));
                                targetEntry(m)=[];
                                m = numel(targetEntry)+1;
                            end
                            m=m-1;
                        end
                    end

                    if numel(targetEntry) ~= numel(targetExit)
                        matRad_cfg.dispError('Inconsistency during ray tracing. Please check correct assignment and overlap priorities of structure types OAR & TARGET.');
                    end

                    beam.ray(j).energy = [];
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
                            raShi.sourceRashiDistance = round(min(ctEntryPoint) - 2*rangeShifterEqD,-1); %place a little away from entry, round to cms to reduce number of unique settings

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

                    beam.ray(j).numParticlesPerMU = numParticlesPerMU;
                    beam.ray(j).minMU = minMU;
                    beam.ray(j).maxMU = maxMU;

                else % target not hit
                    beam.ray(j)               = [];
                    beam.numOfBixelsPerRay(j) = [];
                end
            end
            beam.numOfRays = numel(beam.ray);
        end

        function  beam = finalizeBeam(this,beam)

            % get minimum energy per field
            minEnergy = min([beam.ray.energy]);
            maxEnergy = max([beam.ray.energy]);

            % get corresponding peak position
            minPeakPos  = this.machine.data(minEnergy == this.availableEnergies).peakPos;
            maxPeakPos  = this.machine.data(maxEnergy == this.availableEnergies).peakPos;

            % find set of energyies with adequate spacing
            tolerance              = beam.longitudinalSpotSpacing/10;

            useEnergyBool = this.availablePeakPos >= minPeakPos & this.availablePeakPos <= maxPeakPos;

            ixCurr = find(useEnergyBool,1,'first');
            ixRun  = ixCurr + 1;
            ixEnd  = find(useEnergyBool,1,'last');

            while ixRun <= ixEnd
                if abs(this.availablePeakPos(ixRun)-this.availablePeakPos(ixCurr)) < ...
                        this.longitudinalSpotSpacing - tolerance
                    useEnergyBool(ixRun) = 0;
                else
                    ixCurr = ixRun;
                end
                ixRun = ixRun + 1;
            end

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
            
            %Check additional base data
            available = available && all(isfield(machine.data,{'peakPos','offset'}));

            if ~available
                return;
            else
                available = false;
                msg = [];
            end
    
            %checkBasic
            try    
                %check modality
                checkModality = any(strcmp(matRad_StfGeneratorParticleIMPT.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorParticleIMPT.possibleRadiationModes, pln.radiationMode));
                
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
