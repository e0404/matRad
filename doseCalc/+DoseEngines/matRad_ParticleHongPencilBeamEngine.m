classdef matRad_ParticleHongPencilBeamEngine < DoseEngines.matRad_ParticlePencilBeamEngineAbstract
% matRad_ParticlePencilBeamEngineAbstractGaussian: 
%   Implements an engine for particle based dose calculation 
%   For detailed information see superclass matRad_DoseEngine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the
% help edit

% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
           possibleRadiationModes = {'protons', 'carbon'}
           name = 'Particle Pencil-Beam';
    end
       
    methods 
        
        function this = matRad_ParticleHongPencilBeamEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct
             
            this = this@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(pln);
        end
        
        function dij = calcDose(this,ct,cst,stf)
            % matRad particle dose calculation wrapper
            % can be automaticly called through matRad_calcDose or
            % matRad_calcParticleDose
            %
            % call
            %   dij = this.calcDose(ct,cst,stf,pln)
            %
            % input
            %   ct:             ct cube
            %   cst:            matRad cst struct
            %   stf:            matRad steering information struct
            %
            % output
            %   dij:            matRad dij struct
            %
            % References
            %   [1] http://iopscience.iop.org/0031-9155/41/8/005
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2015 the matRad development team. 
            % 
            % This file is part of the matRad project. It is subject to the license 
            % terms in the LICENSE file found in the top-level directory of this 
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
            % of the matRad project, including this file, may be copied, modified, 
            % propagated, or distributed except according to the terms contained in the 
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            matRad_cfg =  MatRad_Config.instance();
            
            % init dose calc
            [dij,ct,cst,stf] = this.calcDoseInit(ct,cst,stf);
                        
            %Progress Counter
            counter = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:length(stf) % loop over all beams

                % init beam
                dij = this.calcDoseInitBeam(dij,ct,cst,stf,i);                    

                for j = 1:stf(i).numOfRays % loop over all rays

                    if ~isempty(stf(i).ray(j).energy)
                        
                        currRay = this.computeRayGeometry(stf(i).ray(j),dij);

                        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                            counter = counter + 1;
                            this.bixelsPerBeam = this.bixelsPerBeam + 1;

                            this.progressUpdate(counter/dij.totalNumOfBixels);

                            % remember beam and bixel number
                            if ~this.calcDoseDirect
                               dij.beamNum(counter)  = i;
                               dij.rayNum(counter)   = j;
                               dij.bixelNum(counter) = k;

                               % extract MU data if present (checks for downwards compatability)
                                minMU = 0;
                                if isfield(currRay,'minMU')
                                    minMU = currRay.minMU(k);
                                end
    
                                maxMU = Inf;
                                if isfield(currRay,'maxMU')
                                    maxMU = currRay.maxMU(k);
                                end
    
                                numParticlesPerMU = 1e6;
                                if isfield(currRay,'numParticlesPerMU')
                                    numParticlesPerMU = currRay.numParticlesPerMU(k);
                                end
    
                                dij.minMU(counter,1) = minMU;
                                dij.maxMU(counter,1) = maxMU;
                                dij.numParticlesPerMU(counter,1) = numParticlesPerMU;
                            end


                            % find energy index in base data
                            energyIx = find(this.round2(currRay.energy(k),4) == this.round2([this.machine.data.energy],4));

                            % create offset vector to account for additional offsets modelled in the base data and a potential
                            % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                            % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow
                            % interpolations in this.calcParticleDoseBixel() and avoid extrapolations.
                            offsetRadDepth = this.machine.data(energyIx).offset - currRay.rangeShifter(k).eqThickness;

                            % find depth depended lateral cut off
                            if this.dosimetricLateralCutOff == 1
                                currIx = currRay.radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth;
                            elseif this.dosimetricLateralCutOff < 1 && this.dosimetricLateralCutOff > 0
                                % perform rough 2D clipping
                                currIx = currRay.radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                    currRay.radialDist_sq <= max(this.machine.data(energyIx).LatCutOff.CutOff.^2);

                                %currIx = currRay.radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                %    currRay.latDistsX.^2 + currRay.latDistsZ.^2 <= max(this.machine.data(energyIx).LatCutOff.CutOff.^2);

                                % peform fine 2D clipping
                                if length(this.machine.data(energyIx).LatCutOff.CutOff) > 1
                                    currIx(currIx) = matRad_interp1((this.machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                        (this.machine.data(energyIx).LatCutOff.CutOff.^2)', currRay.radDepths(currIx)) >= currRay.radialDist_sq(currIx);
                                end
                            else
                                matRad_cfg.dispError('Cutoff must be a value between 0 and 1!')
                            end

                            % empty bixels may happen during recalculation of error
                            % scenarios -> skip to next bixel
                            if ~any(currIx)
                                continue;
                            end

                            % adjust radDepth according to range shifter
                            currRadDepths = currRay.radDepths(currIx) + currRay.rangeShifter(k).eqThickness;

                            % select correct initial focus sigma squared
                            currSigmaIni_sq = currRay.sigmaIni(k)^2;

                            % consider range shifter for protons if applicable
                            if currRay.rangeShifter(k).eqThickness > 0 && strcmp(stf(i).radiationMode,'protons')

                                % compute!
                                sigmaRashi = matRad_calcSigmaRashi(this.machine.data(energyIx).energy, ...
                                    currRay.rangeShifter(k), ...
                                    currRay.SSD);

                                % add to initial sigma in quadrature
                                currSigmaIni_sq = currSigmaIni_sq +  sigmaRashi^2;

                            end
                            
                            % calculate particle dose for bixel k on ray j of beam i
                            bixelDose = this.calcParticleDoseBixel(...
                                currRadDepths, ...
                                currRay.radialDist_sq(currIx), ...
                                currSigmaIni_sq, ...
                                this.machine.data(energyIx));

                            % dij sampling is not implemented for
                            % particles If we decied to implement it,
                            % it should be implemented as interface in
                            % the superclass as well

                            %{
                                if this.enableDijSampling 
                                    [currIx,bixelDose] = this.dijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);
                                end
                            %}

                            % Save dose for every bixel in cell array
                            this.tmpMatrixContainers.physicalDose{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(currRay.ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);

                            if isfield(dij,'mLETDose')
                                % calculate particle LET for bixel k on ray j of beam i
                                depths = this.machine.data(energyIx).depths + this.machine.data(energyIx).offset;
                                bixelLET = matRad_interp1(depths,this.machine.data(energyIx).LET,currRadDepths);

                                % Save LET for every bixel in cell array
                                this.tmpMatrixContainers.mLETDose{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(currRay.ix(currIx)),1,bixelLET.*bixelDose,dij.doseGrid.numOfVoxels,1);
                            end

                            if this.calcBioDose
                                % calculate alpha and beta values for bixel k on ray j of                  
                                [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                                    currRadDepths,...
                                    currRay.vTissueIndex_j(currIx,:),...
                                    this.machine.data(energyIx));

                                this.tmpMatrixContainers.mAlphaDose{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(currRay.ix(currIx)),1,bixelAlpha.*bixelDose,dij.doseGrid.numOfVoxels,1);
                                this.tmpMatrixContainers.mSqrtBetaDose{mod(counter-1,this.numOfBixelsContainer)+1,1}  = sparse(this.VdoseGrid(currRay.ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.doseGrid.numOfVoxels,1);
                            end

                            %  fill the dij struct each time a
                            %  bixelContainer is calculated and at the end
                            %  of the dose calculation
                            dij = this.fillDij(dij,stf,i,j,k,counter);
                        end
                    end
                end
            end

            %Finalize dose calculation
            dij = this.calcDoseFinalize(ct,cst,stf,dij);
        end

    end

    methods (Access = protected)
        function ray = computeRayGeometry(this,ray,dij)
            ray = computeRayGeometry@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(this,ray,dij);
            ray.currRadialDist_sq = ray.latDistsX.^2 + ray.latDistsZ.^2;
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)   
            % see superclass for information
            
            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');

                %check modality
                checkModality = any(strcmp(DoseEngines.matRad_ParticleHongPencilBeamEngine.possibleRadiationModes, machine.meta.radiationMode));
                
                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            checkMeta = all(isfield(machine.meta,{'SAD','BAMStoIsoDist','LUT_bxWidthminFWHM','dataType'}));

            dataType = machine.meta.dataType;
            if strcmp(dataType,'singleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','sigma','offset','initFocus'}));
            elseif strcmp(dataType,'doubleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','weight','sigma1','sigma2','offset','initFocus'}));
            else
                checkData = false;
            end
            
            available = checkMeta && checkData;
        end
    end
end

