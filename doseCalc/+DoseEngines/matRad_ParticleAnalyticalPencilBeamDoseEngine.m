classdef matRad_ParticleAnalyticalPencilBeamDoseEngine < DoseEngines.matRad_AnalyticalPencilBeamEngine
    % matRad_ParticleDoseEngine: 
    %   Implements an engine for particle based dose calculation 
    %   For detailed information see superclass matRad_DoseEngine
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2015 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the
    % help edit

    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Constant)
           possibleRadiationModes = {'protons', 'carbon'}
           name = 'pencil beam particle';
    end
    
    properties (SetAccess = private, GetAccess = public)  
        
        letDoseTmpContainer; % temporary dose LET container
        alphaDoseTmpContainer; % temporary dose alpha dose container
        betaDoseTmpContainer; % temporary dose beta dose container
        
        calcLET = false; % Boolean which defines if LET should be calculated
        calcBioDose = false; % Boolean which defines if calculation should account for bio optimization
        
        fineSampling; % Boolean switch if using fineSampling  
        fineSamplingN; % number of subsample beams shells, see matRad_calcWeights
        fineSamplingSigmaSub; % gaussian of the sub-beams , see matRad_calcWeights
        fineSamplingMethod; % method used for fine sampling
        
    end
         
    methods 
        
        function obj = matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   ct:                         matRad ct struct
            %   stf:                        matRad steering information struct
            %   pln:                        matRad plan meta information struct
            %   cst:                        matRad cst struct
             
            obj = obj@DoseEngines.matRad_AnalyticalPencilBeamEngine();
            
            if exist('pln','var')
                % check if bio optimization is needed and set the
                % coresponding boolean accordingly
                 if (isfield(pln,'propOpt')&& isfield(pln.propOpt,'bioOptimization')&& ...
                    (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') ||... 
                    isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) && ... 
                    strcmp(pln.radiationMode,'carbon'))
                    obj.calcBioDose = true;
                 end
                 
                 if isfield(pln,'propDoseCalc') && ...
                    isfield(pln.propDoseCalc,'calcLET') && ...
                    pln.propDoseCalc.calcLET
                    obj.calcLET = true;
                 end
                    
            end
        end
        
        function dij = calcDose(obj,ct,stf,pln,cst)
            % matRad particle dose calculation wrapper
            % can be automaticly called through matRad_calcDose or
            % matRad_calcParticleDose
            %
            % call
            %   dij = obj.calcDose(ct,stf,pln,cst)
            %
            % input
            %   ct:             ct cube
            %   stf:            matRad steering information struct
            %   pln:            matRad plan meta information struct
            %   cst:            matRad cst struct
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
            [ct,stf,pln,dij] = obj.calcDoseInit(ct,stf,pln,cst);

            % initialize waitbar
            figureWait = waitbar(0,'calculate dose influence matrix for particles...');
            % prevent closure of waitbar and show busy state
            set(figureWait,'pointer','watch');

            % helper function for energy selection
            round2 = @(a,b)round(a*10^b)/10^b;
            
            % allocate alpha and beta dose container and sparse matrices in the dij struct,
            % for more informations see corresponding method
            dij = obj.allocateBioDoseContainer(dij,pln);
            
            % allocate LET containner and let sparse matrix in dij struct
            if obj.calcLET
                dij = obj.allocateLETContainer(dij,pln);
            end

            % generates tissue class matrix for biological optimization
            if obj.calcBioDose

                if   isfield(obj.machine.data,'alphaX') && isfield(obj.machine.data,'betaX')

                    matRad_cfg.dispInfo('matRad: loading biological base data... ');
                    vTissueIndex = zeros(size(obj.VdoseGrid,1),1);
                    dij.ax       = zeros(dij.doseGrid.numOfVoxels,1);
                    dij.bx       = zeros(dij.doseGrid.numOfVoxels,1);

                    cst = matRad_setOverlapPriorities(cst);

                    % resizing cst to dose cube resolution 
                    cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                                                     dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
                    % retrieve photon LQM parameter for the current dose grid voxels
                    [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,obj.VdoseGrid);

                    for i = 1:size(cst,1)

                        % check if cst is compatiable 
                        if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX') 

                            % check if base data contains alphaX and betaX
                            IdxTissue = find(ismember(obj.machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                                             ismember(obj.machine.data(1).betaX,cst{i,5}.betaX));

                            % check consitency of biological baseData and cst settings
                            if ~isempty(IdxTissue)
                                isInVdoseGrid = ismember(obj.VdoseGrid,cst{i,4}{1});
                                vTissueIndex(isInVdoseGrid) = IdxTissue;
                            else
                                matRad_cfg.dispError('biological base data and cst inconsistent\n');
                            end

                        else
                                vTissueIndex(row) = 1;
                                matRad_cfg.dispInfo(['matRad: tissue type of ' cst{i,2} ' was set to 1\n']);          
                        end
                    end
                    matRad_cfg.dispInfo('done.\n');

                else

                    matRad_cfg.dispError('base data is incomplement - alphaX and/or betaX is missing');

                end

            % issue warning if biological optimization not possible
            elseif sum(strcmp(pln.propOpt.bioOptimization,{'LEMIV_effect','LEMIV_RBExD'}))>0 && ~strcmp(pln.radiationMode,'carbon') ||...
                   ~strcmp(pln.radiationMode,'protons') && strcmp(pln.propOpt.bioOptimization,'const_RBExD')
                warndlg([pln.propOpt.bioOptimization ' optimization not possible with ' pln.radiationMode '- physical optimization is carried out instead.']);
                pln.propOpt.bioOptimization = 'none';      
            end

            % lateral cutoff for raytracing and geo calculations
            obj.effectiveLateralCutoff = matRad_cfg.propDoseCalc.defaultGeometricCutOff;

            matRad_cfg.dispInfo('matRad: Particle dose calculation...\n');
            counter = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:length(stf) % loop over all beams

                % init beam
                dij = obj.calcDoseInitBeam(ct,stf,dij,i);     

                % Determine lateral cutoff
                matRad_cfg.dispInfo('matRad: calculate lateral cutoff...');
                cutOffLevel = matRad_cfg.propDoseCalc.defaultLateralCutOff;
                visBoolLateralCutOff = 0;
                obj.calcLateralParticleCutOff(cutOffLevel,stf(i),visBoolLateralCutOff);
                matRad_cfg.dispInfo('done.\n');    

                for j = 1:stf(i).numOfRays % loop over all rays

                    if ~isempty(stf(i).ray(j).energy)

                        % find index of maximum used energy (round to keV for numerical
                        % reasons
                        energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([obj.machine.data.energy],4);

                        maxLateralCutoffDoseCalc = max(obj.machine.data(energyIx).LatCutOff.CutOff);

                        % calculate initial sigma for all bixel on current ray
                        sigmaIniRay = matRad_calcSigmaIni(obj.machine.data,stf(i).ray(j),stf(i).ray(j).SSD);

                        if strcmp(obj.pbCalcMode, 'fineSampling')
                            % Ray tracing for beam i and ray j
                            [ix,~,~,~,latDistsX,latDistsZ] = obj.calcGeoDists(obj.rot_coordsVdoseGrid, ...
                                                                 stf(i).sourcePoint_bev, ...
                                                                 stf(i).ray(j).targetPoint_bev, ...
                                                                 obj.machine.meta.SAD, ...
                                                                 find(~isnan(obj.radDepthVdoseGrid{1})), ...
                                                                 maxLateralCutoffDoseCalc);

                            % Given the initial sigmas of the sampling ray, this
                            % function provides the weights for the sub-pencil beams,
                            % their positions and their sigma used for dose calculation
                            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
                                if (obj.fineSamplingSigmaSub < sigmaIniRay(k)) && (obj.fineSamplingSigmaSub > 0)
                                    [finalWeight(:,k), sigmaSub(:,k), posX(:,k), posZ(:,k), numOfSub(:,k)] = ...
                                              matRad_calcWeights(sigmaIniRay(k), obj.fineSamplingMethod, obj.fineSamplingN, obj.fineSamplingSigmaSub);
                                else
                                    if (obj.fineSamplingSigmaSub < 0)
                                        matRad_cfg.dispError('Chosen fine sampling sigma cannot be negative!');
                                    elseif (obj.fineSamplingSigmaSub > sigmaIniRay(k))
                                        matRad_cfg.dispError('Chosen fine sampling sigma is too high for defined plan!');
                                    end                          
                                end
                            end
                        else
                            % Ray tracing for beam i and ray j
                            [ix,currRadialDist_sq,~,~,~,~] = obj.calcGeoDists(obj.rot_coordsVdoseGrid, ...
                                                                 stf(i).sourcePoint_bev, ...
                                                                 stf(i).ray(j).targetPoint_bev, ...
                                                                 obj.machine.meta.SAD, ...
                                                                 find(~isnan(obj.radDepthVdoseGrid{1})), ...
                                                                 maxLateralCutoffDoseCalc);

                            radDepths = obj.radDepthVdoseGrid{1}(ix); 
                        end

                        % just use tissue classes of voxels found by ray tracer
                        if obj.calcBioDose
                                vTissueIndex_j = vTissueIndex(ix,:);
                        end



                        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                            counter = counter + 1;
                            obj.bixelsPerBeam = obj.bixelsPerBeam + 1;

                            % Display progress and update text only 200 times
                            if mod(obj.bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                                    matRad_progress(obj.bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                                    floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                            end

                            % update waitbar only 100 times if it is not closed
                            if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                                waitbar(counter/dij.totalNumOfBixels,figureWait);
                            end

                            % remember beam and bixel number
                            if ~obj.calcDoseDirect
                               dij.beamNum(counter)  = i;
                               dij.rayNum(counter)   = j;
                               dij.bixelNum(counter) = k;
                            end

                            % find energy index in base data
                            energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([obj.machine.data.energy],4));


                                if strcmp(obj.pbCalcMode, 'fineSampling')

                                    % calculate projected coordinates for fine sampling of
                                    % each beamlet
                                    projCoords = matRad_projectOnComponents(obj.VdoseGrid(ix), size(obj.radDepthsMat{1}), stf(i).sourcePoint_bev,...
                                                    stf(i).ray(j).targetPoint_bev, stf(i).isoCenter,...
                                                    [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z],...
                                                    -posX(:,k), -posZ(:,k), obj.rotMat_system_T);

                                    % interpolate radiological depths at projected
                                    % coordinates
                                    radDepths = interp3(obj.radDepthsMat{1},projCoords(:,1,:)./dij.doseGrid.resolution.x,...
                                        projCoords(:,2,:)./dij.doseGrid.resolution.y,projCoords(:,3,:)./dij.doseGrid.resolution.z,'nearest');                       

                                    % compute radial distances relative to pencil beam
                                    % component
                                    currRadialDist_sq = reshape(bsxfun(@plus,latDistsX,posX(:,k)'),[],1,numOfSub(k)).^2 + reshape(bsxfun(@plus,latDistsZ,posZ(:,k)'),[],1,numOfSub(k)).^2;
                                end

                                % create offset vector to account for additional offsets modelled in the base data and a potential 
                                % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                                % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow  
                                % interpolations in obj.calcParticleDoseBixel() and avoid extrapolations.
                                offsetRadDepth = obj.machine.data(energyIx).offset - stf(i).ray(j).rangeShifter(k).eqThickness;                               
                                
                                % find depth depended lateral cut off
                                if cutOffLevel >= 1
                                    currIx = radDepths <= obj.machine.data(energyIx).depths(end) + offsetRadDepth;
                                elseif cutOffLevel < 1 && cutOffLevel > 0
                                    % perform rough 2D clipping
                                    currIx = radDepths <= obj.machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                         currRadialDist_sq <= max(obj.machine.data(energyIx).LatCutOff.CutOff.^2);

                                    % peform fine 2D clipping  
                                    if length(obj.machine.data(energyIx).LatCutOff.CutOff) > 1
                                        currIx(currIx) = matRad_interp1((obj.machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                            (obj.machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= currRadialDist_sq(currIx);
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
                                    sigmaRashi = matRad_calcSigmaRashi(obj.machine.data(energyIx).energy, ...
                                                                       stf(i).ray(j).rangeShifter(k), ...
                                                                       stf(i).ray(j).SSD);

                                    % add to initial sigma in quadrature
                                    sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

                                end

                            if strcmp(obj.pbCalcMode, 'fineSampling')
                                % initialise empty dose array
                                totalDose = zeros(size(currIx,1),1);

                                if isfield(dij,'mLETDose')
                                    % calculate particle LET for bixel k on ray j of beam i
                                    depths = obj.machine.data(energyIx).depths + obj.machine.data(energyIx).offset; 
                                    totalLET = zeros(size(currIx,1),1);
                                end

                                % run over components
                                for c = 1:numOfSub
                                    tmpDose = zeros(size(currIx,1),1);
                                    bixelDose = finalWeight(c,k).*obj.calcParticleDoseBixel(...
                                            radDepths(currIx(:,:,c),1,c), ...
                                            currRadialDist_sq(currIx(:,:,c),:,c), ...
                                            sigmaSub(k)^2, ...
                                            obj.machine.data(energyIx));

                                    tmpDose(currIx(:,:,c)) = bixelDose;
                                    totalDose = totalDose + tmpDose;

                                    if isfield(dij,'mLETDose') 
                                        tmpLET = zeros(size(currIx,1),1);
                                        tmpLET(currIx(:,:,c)) = matRad_interp1(depths,obj.machine.data(energyIx).LET,radDepths(currIx(:,:,c),1,c));    
                                        totalLET = totalLET + tmpLET;
                                    end
                                end

                                obj.doseTmpContainer{mod(counter-1,obj.numOfBixelsContainer)+1,1} = sparse(obj.VdoseGrid(ix),1,totalDose,dij.doseGrid.numOfVoxels,1);
                                if isfield(dij,'mLETDose') 
                                    obj.letDoseTmpContainer{mod(counter-1,obj.numOfBixelsContainer)+1,1} = sparse(obj.VdoseGrid(ix),1,totalDose.*totalLET,dij.doseGrid.numOfVoxels,1);
                                end                    
                            else
                                % calculate particle dose for bixel k on ray j of beam i
                                bixelDose = obj.calcParticleDoseBixel(...
                                    currRadDepths, ...
                                    currRadialDist_sq(currIx), ...
                                    sigmaIni_sq, ...
                                    obj.machine.data(energyIx));                 

                                % dij sampling is exluded for particles until we investigated the influence of voxel sampling for particles
                                %relDoseThreshold   =  0.02;   % sample dose values beyond the relative dose
                                %Type               = 'dose';
                                %[currIx,bixelDose] = matRad_DijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);

                                % Save dose for every bixel in cell array
                                obj.doseTmpContainer{mod(counter-1,obj.numOfBixelsContainer)+1,1} = sparse(obj.VdoseGrid(ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);

                                if isfield(dij,'mLETDose')
                                  % calculate particle LET for bixel k on ray j of beam i
                                  depths = obj.machine.data(energyIx).depths + obj.machine.data(energyIx).offset; 
                                  bixelLET = matRad_interp1(depths,obj.machine.data(energyIx).LET,currRadDepths); 

                                  % Save LET for every bixel in cell array
                                  obj.letDoseTmpContainer{mod(counter-1,obj.numOfBixelsContainer)+1,1} = sparse(obj.VdoseGrid(ix(currIx)),1,bixelLET.*bixelDose,dij.doseGrid.numOfVoxels,1);
                                end
                            end



                            if obj.calcBioDose
                                % calculate alpha and beta values for bixel k on ray j of                  
                                [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                                    currRadDepths,...
                                    vTissueIndex_j(currIx,:),...
                                    obj.machine.data(energyIx));

                                obj.alphaDoseTmpContainer{mod(counter-1,obj.numOfBixelsContainer)+1,1} = sparse(obj.VdoseGrid(ix(currIx)),1,bixelAlpha.*bixelDose,dij.doseGrid.numOfVoxels,1);
                                obj.betaDoseTmpContainer{mod(counter-1,obj.numOfBixelsContainer)+1,1}  = sparse(obj.VdoseGrid(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.doseGrid.numOfVoxels,1);
                            end

                            %  fill the dij struct each time a
                            %  bixelContainer is calculated and at the end
                            %  of the dose calculation
                            if mod(counter,obj.numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels                      
                                if obj.calcDoseDirect
                                    dij = obj.fillDijDirect(dij,stf,pln,i,j,k);
                                else
                                    dij = obj.fillDij(dij,stf,pln,counter);
                                end

                            end

                        end

                    end

                end
            end

            %Close Waitbar
            if ishandle(figureWait)
                delete(figureWait);
            end
        end

    end
    
    methods (Access = protected)
        
        function [ct,stf,pln,dij] = calcDoseInit(obj,ct,stf,pln,cst)
            % Extended version of the calcDoseInit method of
            % @matRad_DoseEngine method. See superclass for more information 

            matRad_cfg =  MatRad_Config.instance();
                        
            % assign analytical mode
            if isfield(pln.propDoseCalc,'fineSampling') && strcmp(pln.radiationMode, 'protons')
                obj.pbCalcMode = 'fineSampling';
                defaultFineSampling = matRad_cfg.propDoseCalc.defaultFineSamplingProperties;    
                if isfield(pln.propDoseCalc.fineSampling,'N')
                    obj.fineSamplingN = pln.propDoseCalc.fineSampling.N;
                else
                    obj.fineSamplingN = defaultFineSampling.N;
                end
                if isfield(pln.propDoseCalc.fineSampling,'sigmaSub')    
                    obj.fineSamplingSigmaSub = pln.propDoseCalc.fineSampling.sigmaSub;
                else
                    obj.fineSamplingSigmaSub = defaultFineSampling.sigmaSub;
                end
                if isfield(pln.propDoseCalc.fineSampling,'method')    
                    obj.fineSamplingMethod = pln.propDoseCalc.fineSampling.method;
                else
                    obj.fineSamplingMethod = defaultFineSampling.method;
                end
            else
                obj.pbCalcMode = 'standard';
            end

            [ct,stf,pln,dij] = calcDoseInit@DoseEngines.matRad_DoseEngine(obj,ct,stf,pln,cst);
            
        end
        
        function dij = allocateBioDoseContainer(obj,dij,pln)
        % allocate space for container used in bio optimization

            % get instance of matRad Config for displaying info
            matRad_cfg = MatRad_Config.instance();
            
            if obj.calcBioDose

                    obj.alphaDoseTmpContainer = cell(obj.numOfBixelsContainer,dij.numOfScenarios);
                    obj.betaDoseTmpContainer  = cell(obj.numOfBixelsContainer,dij.numOfScenarios);
                    for i = 1:dij.numOfScenarios
                        dij.mAlphaDose{i}    = spalloc(dij.doseGrid.numOfVoxels,obj.numOfColumnsDij,1);
                        dij.mSqrtBetaDose{i} = spalloc(dij.doseGrid.numOfVoxels,obj.numOfColumnsDij,1);
                    end

            elseif isequal(pln.propOpt.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
                        dij.RBE = 1.1;
                        matRad_cfg.dispInfo('matRad: Using a constant RBE of %g\n',dij.RBE);   
            end
        end
        
        function dij = allocateLETContainer(obj,dij,pln)
        % allocate space for container used in LET calculation
        
              % get MatLab Config instance for displaying warings  
              matRad_cfg = MatRad_Config.instance();
              if isfield(obj.machine.data,'LET')
                obj.letDoseTmpContainer = cell(obj.numOfBixelsContainer,dij.numOfScenarios);
                % Allocate space for dij.dosexLET sparse matrix
                for i = 1:dij.numOfScenarios
                    dij.mLETDose{i} = spalloc(dij.doseGrid.numOfVoxels,obj.numOfColumnsDij,1);
                end
              else
                matRad_cfg.dispWarning('LET not available in the machine data. LET will not be calculated.');
              end
            
        end
                 
        function dij = fillDij(obj,dij,stf,pln,counter)
        % Sequentially fill the sparse matrix dij from the tmpContainer cell array
        %
        % call
        %   dij = fillDij(obj,dij,stf,pln,counter)
        %
        % input
        %   dij:            matRad dij struct
        %   stf:            matRad steering information struct
        %   pln:            matRad plan meta information struct
        %   cst:            counter for indexing current beam, ray and bixel
        %
        % output
        %   dij:            filled dij struct now holding the pre calculated
        %                   dose influence data
        %
        %   see also fillDijDirect
        
            if ~obj.calcDoseDirect
                
                dij.physicalDose{1}(:,(ceil(counter/obj.numOfBixelsContainer)-1)*obj.numOfBixelsContainer+1:counter) = [obj.doseTmpContainer{1:mod(counter-1,obj.numOfBixelsContainer)+1,1}];

                if isfield(dij,'mLETDose')
                    dij.mLETDose{1}(:,(ceil(counter/obj.numOfBixelsContainer)-1)*obj.numOfBixelsContainer+1:counter) = [obj.letDoseTmpContainer{1:mod(counter-1,obj.numOfBixelsContainer)+1,1}];
                end

                if obj.calcBioDose

                    dij.mAlphaDose{1}(:,(ceil(counter/obj.numOfBixelsContainer)-1)*obj.numOfBixelsContainer+1:counter) = [obj.alphaDoseTmpContainer{1:mod(counter-1,obj.numOfBixelsContainer)+1,1}];
                    dij.mSqrtBetaDose{1}(:,(ceil(counter/obj.numOfBixelsContainer)-1)*obj.numOfBixelsContainer+1:counter) = [obj.betaDoseTmpContainer{1:mod(counter-1,obj.numOfBixelsContainer)+1,1}];
                end
            else
                error([dbstack(1).name ' is not intended for direct dose calculation. For filling the dij inside a direct dose calculation please refer to obj.fillDijDirect.']);
            end    
            
        end
        
        function dij = fillDijDirect(obj,dij,stf,pln,currBeamIdx,currRayIdx,currBixelIdx)
        % fillDijDirect - sequentially fill dij, meant for direct calculation only
        %   Fill the sparse matrix physicalDose inside dij with the
        %   indices given by the direct dose calculation
        %   
        %   see also fillDij.      
            if obj.calcDoseDirect
                if isfield(stf(1).ray(1),'weight') && numel(stf(currBeamIdx).ray(currRayIdx).weight) >= currBixelIdx

                    % score physical dose
                    dij.physicalDose{1}(:,currBeamIdx) = dij.physicalDose{1}(:,currBeamIdx) + stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx) * obj.doseTmpContainer{1,1};
                    
                    % write property for mLETDose
                    if isfield(dij,'mLETDose')
                        dij.mLETDose{1}(:,currBeamIdx) = dij.mLETDose{1}(:,currBeamIdx) + stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx) * obj.letDoseTmpContainer{1,1}; 
                    end

                    if obj.calcBioDose
                        % score alpha and beta matrices
                        dij.mAlphaDose{1}(:,currBeamIdx)    = dij.mAlphaDose{1}(:,currBeamIdx) + stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx) * obj.alphaDoseTmpContainer{1,1};
                        dij.mSqrtBetaDose{1}(:,currBeamIdx) = dij.mSqrtBetaDose{1}(:,currBeamIdx) + stf(currBeamIdx).ray(currRayIdx).weight(currBixelIdx) * obj.betaDoseTmpContainer{1,1};

                    end
                else

                    error(['No weight available for beam ' num2str(currBeamIdx) ', ray ' num2str(currRayIdx) ', bixel ' num2str(currBixelIdx)]);

                end
            else
                error([dbstack(1).name 'not available for not direct dose calculation. Refer to obj.fillDij() for a not direct dose calculation.'])
            end
        end
        
        function dose = calcParticleDoseBixel(obj, radDepths, radialDist_sq, sigmaIni_sq, baseData)
        % matRad visualization of two-dimensional dose distributions 
        % on ct including segmentation
        % 
        % call
        %   dose = obj.calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData)
        %
        % input
        %   radDepths:      radiological depths
        %   radialDist_sq:  squared radial distance in BEV from central ray
        %   sigmaIni_sq:    initial Gaussian sigma^2 of beam at patient surface
        %   baseData:       base data required for particle dose calculation
        %
        % output
        %   dose:   particle dose at specified locations as linear vector
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

        % add potential offset
        depths = baseData.depths + baseData.offset;

        % convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
        conversionFactor = 1.6021766208e-02;

        if ~isfield(baseData,'sigma')

            % interpolate depth dose, sigmas, and weights    
            X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma1 baseData.weight baseData.sigma2],radDepths);

            % set dose for query > tabulated depth dose values to zero
            X(radDepths > max(depths),1) = 0;

            % compute lateral sigmas
            sigmaSq_Narr = X(:,2).^2 + sigmaIni_sq;
            sigmaSq_Bro  = X(:,4).^2 + sigmaIni_sq;

            % calculate lateral profile
            L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
            L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
            L = baseData.LatCutOff.CompFac * ((1-X(:,3)).*L_Narr + X(:,3).*L_Bro);

            dose = X(:,1).*L;
        else

            % interpolate depth dose and sigma
            X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma],radDepths);

            %compute lateral sigma
            sigmaSq = X(:,2).^2 + sigmaIni_sq;

            % calculate dose
            dose = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) .* X(:,1) ./(2*pi*sigmaSq);

         end

        % check if we have valid dose values
        if any(isnan(dose)) || any(dose<0)
           error('Error in particle dose calculation.');
        end 
        end
        
        function calcLateralParticleCutOff(obj,cutOffLevel,stf,visBool)
            % matRad function to calculate a depth dependend lateral cutoff 
            % for each pristine particle beam
            % 
            % call
            %   obj.calcLateralParticleCutOff(cutOffLevel,stf,visBool)
            %
            % input
            %   obj:        current engine object includes machine base data file
            %   cutOffLevel:    cut off level - number between 0 and 1
            %   stf:          	matRad steering information struct
            %   visBool:     	toggle visualization (optional)
            %
            % output
            %   machine:    	changes in the object property machine base data file including an additional field representing the lateral
            %                    cutoff
            %   
            % References
            %   -
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

            if cutOffLevel <= 0.98
               warning('a lateral cut off below 0.98 may result in an inaccurate dose calculation') 
            end

            conversionFactor = 1.6021766208e-02;

            % function handle for calculating depth dose for APM
            sumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
                                          exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

            if (cutOffLevel < 0 || cutOffLevel > 1)
               warning('lateral cutoff is out of range - using default cut off of 0.99') 
               cutOffLevel = 0.99;
            end
            % define some variables needed for the cutoff calculation
            vX = [0 logspace(-1,3,1200)]; % [mm]

            % integration steps
            r_mid          = 0.5*(vX(1:end-1) +  vX(2:end))'; % [mm]
            dr             = (vX(2:end) - vX(1:end-1))';
            radialDist_sq  = r_mid.^2;

            % number of depth points for which a lateral cutoff is determined
            numDepthVal    = 35; 

            % helper function for energy selection
            round2 = @(a,b)round(a*10^b)/10^b;

            % extract SSD for each bixel
            vSSD = ones(1,length([stf.ray(:).energy]));
            cnt = 1;
            for i  = 1:length(stf.ray)
                vSSD(cnt:cnt+numel([stf.ray(i).energy])-1) = stf.ray(i).SSD;
                cnt = cnt + numel(stf.ray(i).energy);
            end

            % setup energy, focus index, sigma look up table - only consider unique rows
            [energySigmaLUT,ixUnique]  = unique([[stf.ray(:).energy]; [stf.ray(:).focusIx] ; vSSD]','rows');
            rangeShifterLUT = [stf.ray(:).rangeShifter];
            rangeShifterLUT = rangeShifterLUT(1,ixUnique);

            % find the largest inital beam width considering focus index, SSD and range shifter for each individual energy
            for i = 1:size(energySigmaLUT,1)

                % find index of maximum used energy (round to keV for numerical reasons
                energyIx = max(round2(energySigmaLUT(i,1),4)) == round2([obj.machine.data.energy],4);

                currFoci = energySigmaLUT(i,2);
                sigmaIni = matRad_interp1(obj.machine.data(energyIx).initFocus.dist(currFoci,:)',...
                                          obj.machine.data(energyIx).initFocus.sigma(currFoci,:)',...
                                          energySigmaLUT(i,3));
                sigmaIni_sq = sigmaIni^2;

                % consider range shifter for protons if applicable
                if  strcmp(obj.machine.meta.radiationMode,'protons') && rangeShifterLUT(i).eqThickness > 0  && ~strcmp(obj.machine.meta.machine,'Generic')

                    %get max range shift
                    sigmaRashi = matRad_calcSigmaRashi(obj.machine.data(energyIx).energy, ...
                                                       rangeShifterLUT(i), ...
                                                       energySigmaLUT(i,3));

                    % add to initial sigma in quadrature
                    sigmaIni_sq = sigmaIni_sq +  sigmaRashi.^2;

                end                          

                energySigmaLUT(i,4) = sigmaIni_sq;

            end

            % find for each individual energy the broadest inital beam width
            uniqueEnergies                = unique(energySigmaLUT(:,1));
            largestSigmaSq4uniqueEnergies = NaN * ones(numel(uniqueEnergies),1);
            ix_Max                        = NaN * ones(numel(uniqueEnergies),1);
            for i = 1:numel(uniqueEnergies)
                [largestSigmaSq4uniqueEnergies(i), ix_Max(i)] = max(energySigmaLUT(uniqueEnergies(i) == energySigmaLUT(:,1),4));
            end

            % get energy indices for looping
            vEnergiesIx = find(ismember([obj.machine.data(:).energy],uniqueEnergies(:,1)));
            cnt         = 0;    

            % loop over all entries in the machine.data struct
            for energyIx = vEnergiesIx

                % set default depth cut off - finite value will be set during first iteration
                depthDoseCutOff = inf;

                % get the current integrated depth dose profile
                if isstruct(obj.machine.data(energyIx).Z)
                    idd_org = sumGauss(obj.machine.data(energyIx).depths,obj.machine.data(energyIx).Z.mean,...
                                               obj.machine.data(energyIx).Z.width.^2,...
                                               obj.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd_org = obj.machine.data(energyIx).Z * conversionFactor;
                end

                [~,peakIxOrg] = max(idd_org); 

                % get indices for which a lateral cutoff should be calculated
                cumIntEnergy = cumtrapz(obj.machine.data(energyIx).depths,idd_org);

                peakTailRelation   = 0.5;
                numDepthValToPeak  = ceil(numDepthVal*peakTailRelation);                                                                          % number of depth values from 0 to peak position
                numDepthValTail    = ceil(numDepthVal*(1-peakTailRelation));                                                                      % number of depth values behind peak position
                energyStepsToPeak  = cumIntEnergy(peakIxOrg)/numDepthValToPeak;
                energyStepsTail    = (cumIntEnergy(end)-cumIntEnergy(peakIxOrg))/numDepthValTail;
                % make sure to include 0, peak position and end position
                vEnergySteps       = unique([0:energyStepsToPeak:cumIntEnergy(peakIxOrg) cumIntEnergy(peakIxOrg) ...
                                             cumIntEnergy(peakIxOrg+1):energyStepsTail:cumIntEnergy(end) cumIntEnergy(end)]);

                [cumIntEnergy,ix] = unique(cumIntEnergy);
                depthValues       = matRad_interp1(cumIntEnergy,obj.machine.data(energyIx).depths(ix),vEnergySteps);

                if isstruct(obj.machine.data(energyIx).Z)
                    idd = sumGauss(depthValues,obj.machine.data(energyIx).Z.mean,...
                                               obj.machine.data(energyIx).Z.width.^2,...
                                               obj.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd  = matRad_interp1(obj.machine.data(energyIx).depths,obj.machine.data(energyIx).Z,depthValues) * conversionFactor; 
                end

                cnt = cnt +1 ;
                % % calculate dose in spot
                baseData                   = obj.machine.data(energyIx);
                baseData.LatCutOff.CompFac = 1;   

                for j = 1:numel(depthValues)

                    % save depth value
                    obj.machine.data(energyIx).LatCutOff.depths(j) = depthValues(j);

                    if cutOffLevel == 1
                        obj.machine.data(energyIx).LatCutOff.CompFac   = 1;
                        obj.machine.data(energyIx).LatCutOff.CutOff(j) = Inf;
                    else

                        % calculate dose
                        dose_r = obj.calcParticleDoseBixel(depthValues(j) + baseData.offset, radialDist_sq, largestSigmaSq4uniqueEnergies(cnt), baseData);

                        cumArea = cumsum(2*pi.*r_mid.*dose_r.*dr);
                        relativeTolerance = 0.5; %in [%]
                        if abs((cumArea(end)./(idd(j)))-1)*100 > relativeTolerance
                            warning('LateralParticleCutOff: shell integration is wrong !')
                        end

                        IX = find(cumArea >= idd(j) * cutOffLevel,1, 'first'); 
                        obj.machine.data(energyIx).LatCutOff.CompFac = cutOffLevel^-1;

                        if isempty(IX)
                            depthDoseCutOff = Inf;
                            warning('LateralParticleCutOff: Couldnt find lateral cut off !')
                        elseif isnumeric(IX)
                            depthDoseCutOff = r_mid(IX);
                        end

                        obj.machine.data(energyIx).LatCutOff.CutOff(j) = depthDoseCutOff;

                    end
                end    
            end    

            %% visualization
            if visBool

                % determine which pencil beam should be plotted
                subIx    = ceil(numel(vEnergiesIx)/2);
                energyIx = vEnergiesIx(subIx);

                baseData       = obj.machine.data(energyIx);
                focusIx        = energySigmaLUT(ix_Max(subIx),2);
                maxSSD         = energySigmaLUT(ix_Max(subIx),3);
                rangeShifter   = rangeShifterLUT(ix_Max(subIx));
                TmpCompFac     = baseData.LatCutOff.CompFac;
                baseData.LatCutOff.CompFac = 1;

                % plot 3D cutoff at one specific depth on a rather sparse grid
                sStep         = 0.5;
                vLatX         = -100 : sStep : 100; % [mm]
                dimX          = numel(vLatX);
                midPos        = round(length(vLatX)/2);
                [X,Y]         = meshgrid(vLatX,vLatX);

                radDepths     = [0:sStep:obj.machine.data(energyIx).depths(end)] + obj.machine.data(energyIx).offset;
                radialDist_sq = (X.^2 + Y.^2);
                radialDist_sq = radialDist_sq(:);
                mDose         = zeros(dimX,dimX,numel(radDepths));
                vDoseInt      = zeros(numel(radDepths),1);

                for kk = 1:numel(radDepths)    

                     % calculate initial focus sigma
                     sigmaIni = matRad_interp1(obj.machine.data(energyIx).initFocus.dist(focusIx,:)', ...
                                               obj.machine.data(energyIx).initFocus.sigma(focusIx,:)',maxSSD);
                     sigmaIni_sq = sigmaIni^2;

                     % consider range shifter for protons if applicable
                     if rangeShifter.eqThickness > 0 && strcmp(pln.radiationMode,'protons')

                          % compute!
                          sigmaRashi = matRad_calcSigmaRashi(obj.machine.data(energyIx).energy,rangeShifter,maxSSD);

                          % add to initial sigma in quadrature
                          sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

                     end

                     mDose(:,:,kk) = reshape(obj.calcParticleDoseBixel(radDepths(kk), radialDist_sq, sigmaIni_sq,baseData),[dimX dimX]);

                     [~,IX]           = min(abs((obj.machine.data(energyIx).LatCutOff.depths + obj.machine.data(energyIx).offset) - radDepths(kk)));
                     TmpCutOff        = obj.machine.data(energyIx).LatCutOff.CutOff(IX);    
                     vXCut            = vX(vX<=TmpCutOff);

                     % integration steps
                     r_mid_Cut        = (0.5*(vXCut(1:end-1) +  vXCut(2:end)))'; % [mm]
                     dr_Cut           = (vXCut(2:end) - vXCut(1:end-1))';
                     radialDist_sqCut = r_mid_Cut.^2;    

                     dose_r_Cut       = obj.calcParticleDoseBixel(radDepths(kk), radialDist_sqCut(:), sigmaIni_sq,baseData);

                     cumAreaCut = cumsum(2*pi.*r_mid_Cut.*dose_r_Cut.*dr_Cut);  

                     if ~isempty(cumAreaCut)
                         vDoseInt(kk) = cumAreaCut(end);
                     end
                end

                % obtain maximum dose
                if isstruct(obj.machine.data(energyIx).Z)
                    idd = sumGauss(depthValues,obj.machine.data(energyIx).Z.mean,...
                                               obj.machine.data(energyIx).Z.width.^2,...
                                               obj.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd  = matRad_interp1(obj.machine.data(energyIx).depths,obj.machine.data(energyIx).Z,depthValues) * conversionFactor; 
                end

                [~,peakixDepth] = max(idd); 
                dosePeakPos = obj.calcParticleDoseBixel(obj.machine.data(energyIx).depths(peakixDepth), 0, sigmaIni_sq, baseData);   

                vLevelsDose = dosePeakPos.*[0.01 0.05 0.1 0.9];
                doseSlice   = squeeze(mDose(midPos,:,:));
                figure,set(gcf,'Color',[1 1 1]);
                subplot(311),h=imagesc(squeeze(mDose(midPos,:,:)));hold on;
                set(h,'AlphaData', .8*double(doseSlice>0));
                contour(doseSlice,vLevelsDose,'LevelListMode','manual','LineWidth',2);hold on

                ax = gca;
                ax.XTickLabelMode = 'manual';
                ax.XTickLabel     = strsplit(num2str(ax.XTick*sStep + obj.machine.data(energyIx).offset),' ')';
                ax.YTickLabelMode = 'manual';
                ax.YTickLabel     = strsplit(num2str(ax.YTick*sStep + obj.machine.data(energyIx).offset),' ')';

                plot(1+(obj.machine.data(energyIx).LatCutOff.depths)*sStep^-1,...
                      obj.machine.data(energyIx).LatCutOff.CutOff * sStep^-1 + midPos,'rx');

                legend({'isodose 1%,5%,10% 90%','calculated cutoff'}) ,colorbar,set(gca,'FontSize',12),xlabel('z [mm]'),ylabel('x [mm]');

                entry = obj.machine.data(energyIx);
                if isstruct(entry.Z)
                   idd = sumGauss(entry.depths,entry.Z.mean,entry.Z.width.^2,entry.Z.weight);
                else
                   idd = obj.machine.data(energyIx).Z;
                end
                subplot(312),plot(obj.machine.data(energyIx).depths,idd*conversionFactor,'k','LineWidth',2),grid on,hold on
                             plot(radDepths - obj.machine.data(energyIx).offset,vDoseInt,'r--','LineWidth',2),hold on,
                             plot(radDepths - obj.machine.data(energyIx).offset,vDoseInt * TmpCompFac,'bx','LineWidth',1),hold on,
                legend({'original IDD',['cut off IDD at ' num2str(cutOffLevel) '%'],'cut off IDD with compensation'},'Location','northwest'),
                xlabel('z [mm]'),ylabel('[MeV cm^2 /(g * primary)]'),set(gca,'FontSize',12)     

                totEnergy        = trapz(obj.machine.data(energyIx).depths,idd*conversionFactor) ;
                totEnergyCutOff  = trapz(radDepths,vDoseInt * TmpCompFac) ;
                relDiff          =  ((totEnergy/totEnergyCutOff)-1)*100;   
                title(['rel diff of integral dose ' num2str(relDiff) '%']);
                baseData.LatCutOff.CompFac = TmpCompFac;

                subplot(313),
                if isfield(obj.machine.data(energyIx),'sigma1')
                    yyaxis left;
                    plot(obj.machine.data(energyIx).LatCutOff.depths,obj.machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on
                    plot(obj.machine.data(energyIx).depths,(obj.machine.data(energyIx).sigma1),':','LineWidth',2),grid on,hold on,ylabel('mm')
                    yyaxis right; 
                    plot(obj.machine.data(energyIx).depths,(obj.machine.data(energyIx).sigma2),'-.','LineWidth',2),grid on,hold on,ylabel('mm')
                    legend({'Cutoff','sigma1','sigma2'});
                else
                    yyaxis left;plot(obj.machine.data(energyIx).LatCutOff.depths,obj.machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on,ylabel('mm')
                    yyaxis right;subplot(313),plot(obj.machine.data(energyIx).depths,obj.machine.data(energyIx).sigma,'LineWidth',2),grid on,hold on
                    legend({'Cutoff','sigma'});ylabel('mm')
                end

                set(gca,'FontSize',12),xlabel('z [mm]'),  ylabel('mm')

                % plot cutoff of different energies
                figure,set(gcf,'Color',[1 1 1]);
                cnt = 1;
                for i = vEnergiesIx
                    plot(obj.machine.data(i).LatCutOff.depths,obj.machine.data(i).LatCutOff.CutOff,'LineWidth',1.5),hold on
                    cellLegend{cnt} = [num2str(obj.machine.data(i).energy) ' MeV'];
                    cnt = cnt + 1;
                end
                grid on, grid minor,xlabel('depth in [mm]'),ylabel('lateral cutoff in [mm]')
                title(['cutoff level = ' num2str(cutOffLevel)]),
                ylim = get(gca,'Ylim');    set(gca,'Ylim',[0 ylim(2)+3]),    legend(cellLegend)
            end
        end

        
    end
    
    methods (Static)
        
        function ret = isAvailable(pln)
            % see superclass for information
            ret = any(strcmp(DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine.possibleRadiationModes, pln.radiationMode));
        end
    end
end

