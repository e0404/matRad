classdef matRad_DoseEngineParticlePB < DoseEngines.matRad_DoseEnginePencilBeamParticle
% matRad_DoseEngineParticlePB: 
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
        
        function this = matRad_DoseEngineParticlePB(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct
             
            this = this@DoseEngines.matRad_DoseEnginePencilBeamParticle(pln);
        end
        
        function dij = calcDose(this,ct,cst,stf,pln)
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
            %   pln:            matRad plan meta information struct
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

            % helper function for energy selection
            round2 = @(a,b)round(a*10^b)/10^b;
                        
            %Progress Counter
            counter = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:length(stf) % loop over all beams

                % init beam
                dij = this.calcDoseInitBeam(dij,ct,cst,stf,i);     

                % Determine lateral cutoff
                this.calcLateralParticleCutOff(this.dosimetricLateralCutOff,stf(i));    

                for j = 1:stf(i).numOfRays % loop over all rays

                    if ~isempty(stf(i).ray(j).energy)

                        % find index of maximum used energy (round to keV for numerical
                        % reasons
                        energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([this.machine.data.energy],4);

                        maxLateralCutoffDoseCalc = max(this.machine.data(energyIx).LatCutOff.CutOff);

                        % calculate initial sigma for all bixel on current ray
                        sigmaIniRay = matRad_calcSigmaIni(this.machine.data,stf(i).ray(j),stf(i).ray(j).SSD);

                        
                        % Ray tracing for beam i and ray j
                        [ix,currRadialDist_sq,~,~,~,~] = this.calcGeoDists(this.rot_coordsVdoseGrid, ...
                            stf(i).sourcePoint_bev, ...
                            stf(i).ray(j).targetPoint_bev, ...
                            this.machine.meta.SAD, ...
                            find(~isnan(this.radDepthVdoseGrid{1})), ...
                            maxLateralCutoffDoseCalc);

                        radDepths = this.radDepthVdoseGrid{1}(ix);

                        % just use tissue classes of voxels found by ray tracer
                        if this.calcBioDose
                                vTissueIndex_j = dij.vTissueIndex(ix,:);
                        end

                        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                            counter = counter + 1;
                            this.bixelsPerBeam = this.bixelsPerBeam + 1;

                            % Display progress and update text only 200 times
                            if mod(this.bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                                    matRad_progress(this.bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                                    floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                            end

                            % update waitbar only 100 times if it is not closed
                            if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(this.hWaitbar)
                                waitbar(counter/dij.totalNumOfBixels,this.hWaitbar);
                            end

                            % remember beam and bixel number
                            if ~this.calcDoseDirect
                               dij.beamNum(counter)  = i;
                               dij.rayNum(counter)   = j;
                               dij.bixelNum(counter) = k;

                               % extract MU data if present (checks for downwards compatability)
                                minMU = 0;
                                if isfield(stf(i).ray(j),'minMU')
                                    minMU = stf(i).ray(j).minMU(k);
                                end
    
                                maxMU = Inf;
                                if isfield(stf(i).ray(j),'maxMU')
                                    maxMU = stf(i).ray(j).maxMU(k);
                                end
    
                                numParticlesPerMU = 1e6;
                                if isfield(stf(i).ray(j),'numParticlesPerMU')
                                    numParticlesPerMU = stf(i).ray(j).numParticlesPerMU(k);
                                end
    
                                dij.minMU(counter,1) = minMU;
                                dij.maxMU(counter,1) = maxMU;
                                dij.numParticlesPerMU(counter,1) = numParticlesPerMU;
                            end


                            % find energy index in base data
                            energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([this.machine.data.energy],4));

                            % create offset vector to account for additional offsets modelled in the base data and a potential
                            % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                            % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow
                            % interpolations in this.calcParticleDoseBixel() and avoid extrapolations.
                            offsetRadDepth = this.machine.data(energyIx).offset - stf(i).ray(j).rangeShifter(k).eqThickness;

                            % find depth depended lateral cut off
                            if this.dosimetricLateralCutOff == 1
                                currIx = radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth;
                            elseif this.dosimetricLateralCutOff < 1 && this.dosimetricLateralCutOff > 0
                                % perform rough 2D clipping
                                currIx = radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                    currRadialDist_sq <= max(this.machine.data(energyIx).LatCutOff.CutOff.^2);

                                % peform fine 2D clipping
                                if length(this.machine.data(energyIx).LatCutOff.CutOff) > 1
                                    currIx(currIx) = matRad_interp1((this.machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                        (this.machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= currRadialDist_sq(currIx);
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
                            currRadDepths = radDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness;

                            % select correct initial focus sigma squared
                            sigmaIni_sq = sigmaIniRay(k)^2;

                            % consider range shifter for protons if applicable
                            if stf(i).ray(j).rangeShifter(k).eqThickness > 0 && strcmp(pln.radiationMode,'protons')

                                % compute!
                                sigmaRashi = matRad_calcSigmaRashi(this.machine.data(energyIx).energy, ...
                                    stf(i).ray(j).rangeShifter(k), ...
                                    stf(i).ray(j).SSD);

                                % add to initial sigma in quadrature
                                sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

                            end
                            
                            % calculate particle dose for bixel k on ray j of beam i
                            bixelDose = this.calcParticleDoseBixel(...
                                currRadDepths, ...
                                currRadialDist_sq(currIx), ...
                                sigmaIni_sq, ...
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
                            this.doseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);

                            if isfield(dij,'mLETDose')
                                % calculate particle LET for bixel k on ray j of beam i
                                depths = this.machine.data(energyIx).depths + this.machine.data(energyIx).offset;
                                bixelLET = matRad_interp1(depths,this.machine.data(energyIx).LET,currRadDepths);

                                % Save LET for every bixel in cell array
                                this.letDoseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(ix(currIx)),1,bixelLET.*bixelDose,dij.doseGrid.numOfVoxels,1);
                            end

                            if this.calcBioDose
                                % calculate alpha and beta values for bixel k on ray j of                  
                                [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                                    currRadDepths,...
                                    vTissueIndex_j(currIx,:),...
                                    this.machine.data(energyIx));

                                this.alphaDoseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(ix(currIx)),1,bixelAlpha.*bixelDose,dij.doseGrid.numOfVoxels,1);
                                this.betaDoseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1}  = sparse(this.VdoseGrid(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.doseGrid.numOfVoxels,1);
                            end

                            %  fill the dij struct each time a
                            %  bixelContainer is calculated and at the end
                            %  of the dose calculation
                            if mod(counter,this.numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels                      
                                if this.calcDoseDirect
                                    dij = this.fillDijDirect(dij,stf,pln,i,j,k);
                                else
                                    dij = this.fillDij(dij,stf,pln,counter);
                                end

                            end

                        end

                    end

                end
            end

            %Close Waitbar
            if ishandle(this.hWaitbar)
                delete(this.hWaitbar);
            end
        end

    end            
end

