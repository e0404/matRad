classdef matRad_StfGeneratorParticleSingleBeamlet < matRad_StfGeneratorParticleRayBixelAbstract
% matRad_StfGeneratorPhotonSingleBeamlet: 
%   Creates a single beamlet for particles
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

    properties 
        energy;
        raShiThickness = 50; %Range shifter to be used if useRangeShifter = true;
    end

    properties (Constant)
        name = 'Particle Single Spot';
        shortName = 'ParticleSingleSpot';
        possibleRadiationModes = {'protons','helium','carbon','VHEE'};
    end  
    
    methods 
        function this = matRad_StfGeneratorParticleSingleBeamlet(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorParticleRayBixelAbstract(pln);
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

        function createPatientGeometry(this)
            % Simple patient geometry with isocenter target
            matRad_cfg = MatRad_Config.instance();

            if isempty(this.isoCenter)
                this.isoCenter = matRad_getIsoCenter(this.cst,this.ct,visBool);
            end

            if ~isequal(size(this.isoCenter),[this.numOfBeams,3]) && ~size(this.isoCenter,1) ~= 1
                matRad_cfg.dispWarning('IsoCenter invalid, creating new one automatically!');
                this.isoCenter = matRad_getIsoCenter(this.cst,this.ct,visBool);
            end
            
            if size(this.isoCenter,1) == 1          
                this.isoCenter = repmat(this.isoCenter,this.numOfBeams,1);
            end            

            %Voxel index of Isocenter
            isoIx = matRad_world2cubeIndex(this.isoCenter,this.ct,true);
            isoIx(isoIx < 1) = 1;
            for i = 1:size(isoIx,1)
                isoIx(i,isoIx(i,:) > this.ct.cubeDim) = this.ct.cubeDim(isoIx(i,:) > this.ct.cubeDim);
            end

            % generate voi cube for targets
            this.voiTarget    = zeros(this.ct.cubeDim);
            for i = 1:size(isoIx,1)
                this.voiTarget(isoIx(i,1),isoIx(i,2),isoIx(i,3)) = 1;
            end

            % Margin info
            if this.addMargin
                adds = unique(perms([1 0 0]),'rows');
                adds = [adds; -adds];
                for i = 1:size(isoIx,1)
                    for p = 1:size(adds,1)
                        padIx = isoIx(i,:) + adds(p,:);
                        if any(padIx < 1) || any(padIx > this.ct.cubeDim)
                            continue;
                        end
                        this.voiTarget(padIx(1),padIx(2),padIx(3)) = 1;
                    end
                end
            end

            V = find(this.voiTarget > 0);

            % throw error message if no target is found
            if isempty(V)
                matRad_cfg.dispError('Could not identify find isocenter target.');
            end           

            % Convert linear indices to 3D voxel coordinates
            this.voxTargetWorldCoords = matRad_cubeIndex2worldCoords(V, this.ct);

            % take only voxels inside patient
            V = [this.cst{:,4}];
            V = unique(vertcat(V{:}));

            % ignore densities outside of contours
            eraseCtDensMask = ones(prod(this.ct.cubeDim), 1);
            eraseCtDensMask(V) = 0;
            for i = 1:this.ct.numOfCtScen
                this.ct.cube{i}(eraseCtDensMask == 1) = 0;
            end
        end

        function beam = setBeamletEnergies(this,beam)
            isoCenterCubeSystem = matRad_world2cubeCoords(beam.isoCenter,this.ct);

            % ray tracing necessary to determine depth of the target
            [alphas,l,rho,d12,~] = matRad_siddonRayTracer(isoCenterCubeSystem, ...
                this.ct.resolution, ...
                beam.sourcePoint, ...
                beam.ray.targetPoint, ...
                [{this.ct.cube{1}} {this.voiTarget}]);

            if isempty(alphas)
                matRad_cfg.dispError('Beam seems to not hit the CT! Check Isocenter placement!');
            end

            ctEntryPoint = alphas(1) * d12;

            if sum(rho{2}) <= 0 && isempty(this.energy)
                availSorted = sort(this.availableEnergies);
                useEnergy = ceil(availSorted(numel(availSorted)/2));
                matRad_cfg.dispWarning('Could not obtain suitable energy based on isoCenter. Taking median energy of %g MeV',useEnergy);
            elseif ~isempty(this.energy)
                [~,ix] = min(abs(this.energy-this.availableEnergies));                
                useEnergy = this.availableEnergies(ix);
            else
                % compute radiological depths
                % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                radDepths = cumsum(l .* rho{1});

                % find target entry & exit
                diff_voi    = diff([rho{2}]);
                if rho{2}(1) > 0
                    targetEntry = 1;
                else
                    targetEntry = radDepths(diff_voi == 1);
                end

                targetExit  = radDepths(diff_voi == -1);

                if numel(targetEntry) ~= numel(targetExit)
                    matRad_cfg.dispError('Inconsistency during ray tracing. Please check correct assignment and overlap priorities of structure types OAR & TARGET.');
                end

                % Save energies in stf struct
                bestPeakPos = mean([targetExit,targetEntry]);
                if this.useRangeShifter
                    bestPeakPos = bestPeakPos + this.raShiThickness;
                end                    

                [~,closest] = min(abs(this.availablePeakPos - bestPeakPos));
                useEnergy = this.availableEnergies(closest);                
            end

            beam.ray.energy = useEnergy;
            % book keeping energy index
            [~, vEnergyIx] = min(abs(beam.ray.energy-this.availableEnergies));

            % get the focus index
            if isfield(this.machine.meta,'LUT_bxWidthminFWHM')
                LUTspotSize = this.machine.meta.LUT_bxWidthminFWHM;
            elseif isfield(this.machine.meta,'LUTspotSize')
                LUTspotSize = this.machine.meta.LUTspotSize;
            else
                LUTspotSize = [];
            end
                    
            if ~isempty(LUTspotSize)
                currentMinimumFWHM = matRad_interp1(LUTspotSize(1,:)',...
                    LUTspotSize(2,:)',...
                    beam.bixelWidth, ...
                    LUTspotSize(2,end));
                
                beam.ray.focusIx = find(this.machine.data(vEnergyIx).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
            else
                beam.ray.focusIx = 1;
            end
            
            % Get Range shifter
            if this.useRangeShifter
            
                %Include range shifter data
                beam.ray.rangeShifter.ID = 1;
                beam.ray.rangeShifter.eqThickness = this.raShiThickness;
            
                %Place range shifter 2 times the range away from isocenter, but
                %at least 10 cm
                sourceRaShi = round(ctEntryPoint - 2*this.raShiThickness,-1); %place a little away from entry, round to cms to reduce number of unique settings;
                beam.ray.rangeShifter.sourceRashiDistance = sourceRaShi;
            else
                beam.ray.rangeShifter.ID = 0;
                beam.ray.rangeShifter.eqThickness = 0;
                beam.ray.rangeShifter.sourceRashiDistance = 0;
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
            [available,msg] = matRad_StfGeneratorParticleRayBixelAbstract.isAvailable(pln,machine);

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
                checkModality = any(strcmp(matRad_StfGeneratorParticleSingleBeamlet.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_StfGeneratorParticleSingleBeamlet.possibleRadiationModes, pln.radiationMode));
                
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
