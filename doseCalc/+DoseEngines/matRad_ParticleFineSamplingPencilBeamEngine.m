classdef matRad_ParticleFineSamplingPencilBeamEngine < DoseEngines.matRad_ParticlePencilBeamEngineAbstract
% matRad_ParticlePencilBeamEngineAbstractFineSampling: 
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
           name = 'Subsampling Particle Pencil-Beam';
    end
    
    properties (SetAccess = public, GetAccess = public)
        fineSampling;                   % Struct with finesampling properties
    end
                 
    methods 
        
        function this = matRad_ParticleFineSamplingPencilBeamEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct
             
            this = this@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(pln);
        end

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(this);

            matRad_cfg = MatRad_Config.instance();
            this.fineSampling = matRad_cfg.propDoseCalc.defaultFineSamplingProperties;
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
                        [ix,~,~,~,latDistsX,latDistsZ] = this.calcGeoDists(this.rot_coordsVdoseGrid, ...
                            stf(i).sourcePoint_bev, ...
                            stf(i).ray(j).targetPoint_bev, ...
                            this.machine.meta.SAD, ...
                            find(~isnan(this.radDepthVdoseGrid{1})), ...
                            maxLateralCutoffDoseCalc);

                        % Given the initial sigmas of the sampling ray, this
                        % function provides the weights for the sub-pencil beams,
                        % their positions and their sigma used for dose calculation
                        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
                            if (this.fineSampling.sigmaSub < sigmaIniRay(k)) && (this.fineSampling.sigmaSub > 0)
                                [finalWeight(:,k), sigmaSub(:,k), posX(:,k), posZ(:,k), numOfSub(:,k)] = ...
                                    this.calcFineSamplingMixture(sigmaIniRay(k));
                            else
                                if (this.fineSamplingSigmaSub < 0)
                                    matRad_cfg.dispError('Chosen fine sampling sigma cannot be negative!');
                                elseif (this.fineSamplingSigmaSub > sigmaIniRay(k))
                                    matRad_cfg.dispError('Chosen fine sampling sigma is too high for defined plan!');
                                end
                            end
                        end
                        

                        % just use tissue classes of voxels found by ray tracer
                        if this.calcBioDose
                                vTissueIndex_j = vTissueIndex(ix,:);
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

                            % calculate projected coordinates for fine sampling of
                            % each beamlet
                            projCoords = matRad_projectOnComponents(this.VdoseGrid(ix), size(this.radDepthCube{1}), stf(i).sourcePoint_bev,...
                                stf(i).ray(j).targetPoint_bev, stf(i).isoCenter,...
                                [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z],...
                                -posX(:,k), -posZ(:,k), this.rotMat_system_T);

                            % interpolate radiological depths at projected
                            % coordinates
                            radDepths = interp3(this.radDepthCube{1},projCoords(:,1,:)./dij.doseGrid.resolution.x,...
                                projCoords(:,2,:)./dij.doseGrid.resolution.y,projCoords(:,3,:)./dij.doseGrid.resolution.z,'nearest');

                            % compute radial distances relative to pencil beam
                            % component
                            currRadialDist_sq = reshape(bsxfun(@plus,latDistsX,posX(:,k)'),[],1,numOfSub(k)).^2 + reshape(bsxfun(@plus,latDistsZ,posZ(:,k)'),[],1,numOfSub(k)).^2;

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

                            % initialise empty dose array
                            totalDose = zeros(size(currIx,1),1);

                            if isfield(dij,'mLETDose')
                                % calculate particle LET for bixel k on ray j of beam i
                                depths = this.machine.data(energyIx).depths + this.machine.data(energyIx).offset;
                                totalLET = zeros(size(currIx,1),1);
                            end

                            % run over components
                            for c = 1:numOfSub(k)
                                tmpDose = zeros(size(currIx,1),1);
                                bixelDose = finalWeight(c,k).*this.calcParticleDoseBixel(...
                                    radDepths(currIx(:,:,c),1,c), ...
                                    currRadialDist_sq(currIx(:,:,c),:,c), ...
                                    sigmaSub(k)^2, ...
                                    this.machine.data(energyIx));

                                tmpDose(currIx(:,:,c)) = bixelDose;
                                totalDose = totalDose + tmpDose;

                                if isfield(dij,'mLETDose')
                                    tmpLET = zeros(size(currIx,1),1);
                                    tmpLET(currIx(:,:,c)) = matRad_interp1(depths,this.machine.data(energyIx).LET,radDepths(currIx(:,:,c),1,c));
                                    totalLET = totalLET + tmpLET;
                                end
                            end

                            this.doseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(ix),1,totalDose,dij.doseGrid.numOfVoxels,1);
                            if isfield(dij,'mLETDose')
                                this.letDoseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(ix),1,totalDose.*totalLET,dij.doseGrid.numOfVoxels,1);
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
    
    methods (Access = protected)

        function dij = calcDoseInitBeam(this,dij,ct,cst,stf,i)
            % Method for initializing the beams for analytical pencil beam
            % dose calculation
            %
            % call
            %   this.calcDoseInitBeam(dij,ct,cst,stf,i)
            %
            % input
            %   dij:                        matRad dij struct
            %   ct:                         matRad ct struct
            %   cst:                        matRad cst struct
            %   stf:                        matRad steering information struct
            %   i:                          index of beam
            %
            % output
            %   dij:                        updated dij struct

            if ~this.keepRadDepthCubes
                this.keepRadDepthCubes = true;
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispInfo('Keeping radiological depth cubes for fine-sampling!');
            end

            dij = calcDoseInitBeam@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(this,dij,ct,cst,stf,i);

        end
                                                    
        function [finalWeight, sigmaBeamlet, posX, posY, numOfSub] = calcFineSamplingMixture(this,sigmaTot)
            % This function creates a Gaussian Mixture Model on a Gaussian
            % for Fine-Sampling
            %
            % call
            %   [finalWeight, sigmaBeamlet, posX, posY, numOfSub] = ...
            %                       this.calcFineSamplingMixture(sigmaTot)
            %
            % input
            %   sigmaTot:       the standard deviation of the lateral spread of the pencil
            %                   beam
            % output
            %   finalWeight:    is the array of the weights of the sub-pencil beams. It
            %                   runs over the same index as posx and posy
            %
            %   posX & posY:    are the positions of the sub-pencil beams, returned as
            %                   meshgrid-matrices if method is 'square' and as vectors
            %                   if method is 'circle'
            %   numOfSub:       number of sub-pencil beams
            %
            % References
            %   [1] https://iopscience.iop.org/article/10.1088/0031-9155/61/1/183
            

            % method:   method of weight calculation
            %           'russo' for integral calculation according to [1]
            %           'fitCircle'   for using square grid weights derived from a fit
            %           'fitSquare'   for using circular grid weights derived from a fit
            method = this.fineSampling.method;
            
            %   N:              if method == russo:
            %                   number of subsample beams shells. Means we have a
            %                   grid of NxN sub beams representing the total beam
            %                   if method == fitCircle or fitSquare:
            %                   number of subsample beams shells. n = 2 means we have two
            %                   lines of points around the central ray. The number of
            %                   sub-beams will be:
            %                   #sb = (2*n +1)^2 for the square;
            %                   #sb = (2^n -1)*6 +1 for the circle
            N = this.fineSampling.N;

            %   sigmaSub:       is the sigma of the gaussian of the sub-beams (only
            %                   needed for russo method)
            %
            sigmaSub = this.fineSampling.sigmaSub;

            if ~strcmp(method, 'russo') && ~strcmp(method, 'fitCircle') && ~strcmp(method, 'fitSquare')

                error('method not supported');

            elseif strcmp(method, 'russo')
                % splitting into N^2 sub beams accoring to Russo et al (2016)

                sigmaHead = sqrt(sigmaTot^2 - sigmaSub^2);

                gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2));

                R = 3.5 * sigmaHead;
                dR = 2 * R / N;

                counter = 1;
                for iX = -(N - 1) / 2 : (N - 1) / 2
                    for iY = -(N - 1) / 2 : (N - 1) / 2
                        finalWeight(counter) = integral2(@(x,y) gauss(sigmaHead, x, y, 0, 0), ...
                            (iX - 0.5) * dR, (iX + 0.5) * dR, ...
                            (iY - 0.5) * dR, (iY + 0.5) * dR);
                        posX(counter) = iX * dR;
                        posY(counter) = iY * dR;
                        sigmaBeamlet(counter) = sigmaSub;

                        counter = counter + 1;
                    end
                end

                finalWeight = finalWeight';
                posX = posX';
                posY = posY';
                sigmaBeamlet = sigmaBeamlet';


                numOfSub = N * N;

            elseif strcmp(method, 'fitCircle') || strcmp(method, 'fitSquare')
                % number of sub beams will be (2*n +1)^2 for the square;
                %                             (2^n -1)*6 +1 for the circle
                if N~=2 && N~=3 && N~=8
                    error('number of shells N not supported');
                end

                % This parameters come from simulations done previously
                % see "Research on the dosimetric accuracy of pencil beam fine sampling
                % for radiation proton dose calculation" by Giuseppe Pezzano (2018)
                if N == 2
                    if strcmp(method,'fitCircle')
                        sigmaBeamlet = 0.8237 .* sigmaTot;
                        radius    = 0.6212 .* sigmaTot;
                        X1(1,:)     = 0.3866 .* sigmaTot.^2;
                        X1(2,:)     = 0.6225 .* sigmaTot;
                    elseif strcmp(method,'fitSquare')
                        sigmaBeamlet = 0.8409 .* sigmaTot;
                        radius    = 0.5519 .* sigmaTot;
                        X1(1,:)     = 0.3099 .* sigmaTot.^2;
                        X1(2,:)     = 0.5556 .* sigmaTot;
                    end
                elseif N == 3
                    if strcmp(method,'fitCircle')
                        sigmaBeamlet = 0.7605 .* sigmaTot;
                        radius    = 0.5000 .* sigmaTot;
                        X1(1,:)     = 0.3006 .* sigmaTot.^2 - 1.3005 .* sigmaTot + 7.3097;
                        X1(2,:)     = 0.6646 .* sigmaTot - 0.0044;
                    elseif strcmp(method,'fitSquare')
                        sigmaBeamlet = 0.8409 .* sigmaTot;
                        radius    = 0.5391 .* sigmaTot + 0.0856;
                        X1(1,:)     = 0.3245 .* sigmaTot.^2 + 0.0001 .* sigmaTot - 0.0004;
                        X1(2,:)     = 0.6290 .* sigmaTot - 0.0403;
                    end
                elseif N == 8 && strcmp(method,'fitCircle')
                    sigmaBeamlet = 0.5 .* sigmaTot;
                    radius    = 0.25 .* sigmaTot;
                    X1(1,:)     = 0.0334 .* sigmaTot.^2 - 4.1061e-06 .* sigmaTot + 1.5047e-06;
                    X1(2,:)     = 0.6 .* sigmaTot + 3.3151e-06;
                else
                    error('number of shells N not supported');
                end

                % setting positions of sub-beams
                if strcmp(method,'fitSquare')
                    numOfSub = (2*N +1)^2;
                    points   = linspace(-radius*N,radius*N,2*N +1);
                    posX     = points'*ones(1,2*N +1);
                    posY     = posX';
                else
                    dim = size(radius,2);
                    numOfSub = (2^N -1)*6 +1;
                    ang  = zeros(1,1);
                    posX = zeros(1,dim);
                    posY = zeros(1,dim);
                    radiusShell = zeros(1,dim);
                    for i = 1:N
                        subsInShell = 6 * 2^(i-1);
                        % this takes the sub-beams index in one shell
                        ang         = cat(2, ang, pi .* linspace(0,2-2/subsInShell, subsInShell));
                        radiusShell = cat(1, radiusShell, ones(subsInShell,1)*(i.*radius));
                    end
                    posX = cat(1, posX, bsxfun(@times,cos(ang(2:end))',radiusShell(2:end,:)));
                    posY = cat(1, posY, bsxfun(@times,sin(ang(2:end))',radiusShell(2:end,:)));
                end

                % compute weights at positions
                sig  = ones(size(posX,1),1)*X1(2,:);
                normSig = ones(size(posX,1),1)*X1(1,:);

                finalWeight = normSig .* (2.*pi.*sig.^2).^(-1) .* exp(-(posX.^2+posY.^2)./(2.*(sig.^2)));
                finalWeight = reshape(finalWeight, numel(finalWeight), 1);
                sigmaBeamlet = repmat(reshape(sigmaBeamlet, numel(sigmaBeamlet), 1), numel(finalWeight),1);
                posX =reshape(posX, numel(posX), 1);
                posY =reshape(posY, numel(posY), 1);

            end

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
                checkModality = any(strcmp(DoseEngines.matRad_DoseEngineParticlePB.possibleRadiationModes, machine.meta.radiationMode));
                
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

