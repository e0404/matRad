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
            this.fineSampling.method = 'fitCircle';
            this.fineSampling.N = 2;
        end
    end
    
    methods (Access = protected)

        function ray = initRay(this,beam,j)           
            ray = initRay@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(this,beam,j);
            
            %We need some more beam information here
            ray.rotMat_system_T = beam.rotMat_system_T;
        end

        % We override this function to get full lateral distances
        function ray = getRayGeometryFromBeam(this,ray,currBeam)
            lateralRayCutOff = this.getLateralDistanceFromDoseCutOffOnRay(ray);

            % Ray tracing for beam i and ray j
            [ray.ix,ray.radialDist_sq,ray.latDists] = this.calcGeoDists(currBeam.rot_coordsVdoseGrid, ...
                ray.sourcePoint_bev, ...
                ray.targetPoint_bev, ...
                ray.SAD, ...
                currBeam.ixRadDepths, ...
                lateralRayCutOff);
            
            ray.radDepths = currBeam.radDepthVdoseGrid{1}(ray.ix);
        end

        function bixel = initBixel(this,currRay,k)
            bixel = initBixel@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(this,currRay,k);
            
            bixel.latDists = currRay.latDists(bixel.subRayIx,:);            

            % Given the initial sigmas of the sampling ray, this
            % function provides the weights for the sub-pencil beams,
            % their positions and their sigma used for dose calculation
            if (this.fineSampling.sigmaSub^2 < bixel.sigmaIniSq) && (this.fineSampling.sigmaSub > 0)
                [bixel.finalSubWeight, bixel.sigmaSub, bixel.subPosX, bixel.subPosZ, bixel.numOfSub] = ...
                    this.calcFineSamplingMixture(sqrt(bixel.sigmaIniSq));
            else
                matRad_cfg = MatRad_Config.instance();
                if (this.fineSampling.sigmaSub < 0)
                    matRad_cfg.dispError('Chosen fine sampling sigma cannot be negative!');
                elseif (this.fineSampling.sigmaSub > sqrt(bixel.sigmaIniSq))
                    matRad_cfg.dispError('Chosen fine sampling sigma is too high for defined plan!');
                end
            end

            % calculate projected coordinates for fine sampling of
            % each beamlet
            projCoords = matRad_projectOnComponents(this.VdoseGrid(bixel.ix), size(this.radDepthCubes{currRay.beamIndex}), currRay.sourcePoint_bev,...
                currRay.targetPoint_bev, currRay.isoCenter,...
                [this.doseGrid.resolution.x this.doseGrid.resolution.y this.doseGrid.resolution.z],...
                -bixel.subPosX, -bixel.subPosZ, currRay.rotMat_system_T);


            % interpolate radiological depths at projected
            % coordinates
            % TODO: we get NaN's here - why? (They come from the
            % radDepthCubes, but I don't know why we interpolate at these
            % positions)
            bixel.radDepths = interp3(this.radDepthCubes{currRay.beamIndex},projCoords(:,1,:)./this.doseGrid.resolution.x,...
                projCoords(:,2,:)./this.doseGrid.resolution.y,projCoords(:,3,:)./this.doseGrid.resolution.z,'nearest',0);


            % compute radial distances relative to pencil beam
            % component
            %bixel.radialDist_sq = reshape(bsxfun(@plus,bixel.latDists(:,1),bixel.subPosX'),[],1,bixel.numOfSub).^2 + reshape(bsxfun(@plus,bixel.latDists(:,2),bixel.subPosZ'),[],1,bixel.numOfSub).^2;
            bixel.radialDist_sq = reshape((bixel.latDists(:,1) + bixel.subPosX').^2 + (bixel.latDists(:,2) + bixel.subPosZ').^2,[],1,bixel.numOfSub);
        end

        function currBeam = initBeam(this,dij,ct,cst,stf,i)
            % Method for initializing the beams for analytical pencil beam
            % dose calculation
            %
            % call
            %   this.initBeam(dij,ct,cst,stf,i)
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

            currBeam = initBeam@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(this,dij,ct,cst,stf,i);
        end
        
        function kernels = interpolateKernelsInDepth(this,bixel)
            kernels = interpolateKernelsInDepth@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(this,bixel);
            kernels = structfun(@(x) reshape(x,size(bixel.radDepths)),kernels,'UniformOutput',false);
        end
        
        function bixel = calcParticleBixel(this,bixel)
            kernel = this.interpolateKernelsInDepth(bixel);
            % initialise empty dose array
            bixel.physicalDose = zeros(size(bixel.ix,1),1);

            if this.calcLET
                bixel.mLETDose = zeros(size(bixel.physicalDose));
            end
            
            %We only have one bixel
            if ~isfield(bixel,'numOfSub')
                bixel.numOfSub = 1;
                bixel.finalSubWeight = 1;
                bixel.sigmaSub = sqrt(bixel.sigmaIniSq);
            end
            
            bixel.sigmaSubSq = bixel.sigmaSub.^2;            
            
            dg = ~isfield(bixel.baseData,'sigma');
      
            if dg
                % compute lateral sigmas
                sigmaSqNarrow = squeeze(kernel.sigma1).^2 + bixel.sigmaSubSq';
                sigmaSqBroad  = squeeze(kernel.sigma2).^2 + bixel.sigmaSubSq';

                % calculate lateral profile
                L_Narr =  exp( -squeeze(bixel.radialDist_sq) ./ (2*sigmaSqNarrow))./(2*pi*sigmaSqNarrow);
                L_Bro  =  exp( -squeeze(bixel.radialDist_sq) ./ (2*sigmaSqBroad ))./(2*pi*sigmaSqBroad );
                L = (1-squeeze(kernel.weight)).*L_Narr + squeeze(kernel.weight).*L_Bro;
            else
                %compute lateral sigma
                sigmaSq = squeeze(kernel.sigma).^2 + bixel.sigmaSubSq';
                L = exp( -squeeze(bixel.radialDist_sq) ./ (2*sigmaSq))./ (2*pi*sigmaSq);
            end
            
            tmpDose = (L .* squeeze(kernel.Z));

            bixel.physicalDose = bixel.baseData.LatCutOff.CompFac*(tmpDose*bixel.finalSubWeight);

            nanIx = isnan(bixel.physicalDose);
            bixel.pyhsicalDose(nanIx) = 0;
            if this.calcLET
                tmpLET = bixel.baseData.LatCutOff.CompFac*((tmpDose .* squeeze(kernel.LET)) *bixel.finalSubWeight);
                bixel.mLETDose(~nanIx) = bixel.mLETDose(~nanIx) + tmpLET(~nanIx);
            end
            
            if this.calcBioDose
                bixel.mAlphaDose = bixel.physicalDose;
                bixel.mSqrtBetaDose = bixel.physicalDose;
                %From matRad_calcLQParameter
                numOfTissueClass = size(bixel.baseData.alpha,2);
                for i = 1:numOfTissueClass
                    mask = bixel.vTissueIndex == i;
                    if any(mask)
                        bixel.mAlphaDose(mask) = bixel.mAlphaDose(mask) .* X.alpha(mask);
                        bixel.mSqrtBetaDose(mask)  = bixel.mSqrtBetaDose(mask) .* X.beta(mask);
                    end
                end
            end
            
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
                        if (iX*dR)^2 + (iY*dR)^2 > R^2
                            continue;
                        end

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
                finalWeight = finalWeight*1/sum(finalWeight);
                posX = posX';
                posY = posY';
                sigmaBeamlet = sigmaBeamlet';


                numOfSub = numel(finalWeight);

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
                checkModality = any(strcmp(DoseEngines.matRad_ParticleFineSamplingPencilBeamEngine.possibleRadiationModes, machine.meta.radiationMode));
                
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
            
            available = checkMeta & checkData;
        end
    end
end

