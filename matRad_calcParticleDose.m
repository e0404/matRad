function dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect)
% matRad particle dose calculation wrapper
% 
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   calcDoseDirect: boolian switch to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
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
matRad_calcDoseInit;

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for particles...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
        && strcmp(pln.radiationMode,'carbon')
   
        alphaDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
        betaDoseTmpContainer  = cell(numOfBixelsContainer,dij.numOfScenarios);
        for i = 1:dij.numOfScenarios
            dij.mAlphaDose{i}    = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
            dij.mSqrtBetaDose{i} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
        end
        
elseif isequal(pln.propOpt.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
            dij.RBE = 1.1;
            matRad_cfg.dispInfo('matRad: Using a constant RBE of %g\n',dij.RBE);   
end

if isfield(pln,'propDoseCalc') && ...
   isfield(pln.propDoseCalc,'calcLET') && ...
   pln.propDoseCalc.calcLET
  if isfield(machine.data,'LET')
    letDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
    % Allocate space for dij.dosexLET sparse matrix
    for i = 1:dij.numOfScenarios
        dij.mLETDose{i} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
    end
  else
    matRad_cfg.dispWarning('LET not available in the machine data. LET will not be calculated.');
  end
end

% generates tissue class matrix for biological optimization
if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
        && strcmp(pln.radiationMode,'carbon')
    
    if   isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
            
        matRad_cfg.dispInfo('matRad: loading biological base data... ');
        vTissueIndex = zeros(size(VdoseGrid,1),1);
        dij.ax       = zeros(dij.doseGrid.numOfVoxels,1);
        dij.bx       = zeros(dij.doseGrid.numOfVoxels,1);

        cst = matRad_setOverlapPriorities(cst);
    
        % resizing cst to dose cube resolution 
        cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                                         dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
        % retrieve photon LQM parameter for the current dose grid voxels
        [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);

        for i = 1:size(cst,1)

            % check if cst is compatiable 
            if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX') 

                % check if base data contains alphaX and betaX
                IdxTissue = find(ismember(machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                                 ismember(machine.data(1).betaX,cst{i,5}.betaX));

                % check consitency of biological baseData and cst settings
                if ~isempty(IdxTissue)
                    isInVdoseGrid = ismember(VdoseGrid,cst{i,4}{1});
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
effectiveLateralCutoff = matRad_cfg.propDoseCalc.defaultGeometricCutOff;

matRad_cfg.dispInfo('matRad: Particle dose calculation...\n');
counter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(stf) % loop over all beams
  
    % init beam
    matRad_calcDoseInitBeam;     
        
    % Determine lateral cutoff
    matRad_cfg.dispInfo('matRad: calculate lateral cutoff...');
    cutOffLevel = matRad_cfg.propDoseCalc.defaultLateralCutOff;
    visBoolLateralCutOff = 0;
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
    matRad_cfg.dispInfo('done.\n');    

    for j = 1:stf(i).numOfRays % loop over all rays

        if ~isempty(stf(i).ray(j).energy)

            % find index of maximum used energy (round to keV for numerical
            % reasons
            energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);

            maxLateralCutoffDoseCalc = max(machine.data(energyIx).LatCutOff.CutOff);
    
            % calculate initial sigma for all bixel on current ray
            sigmaIniRay = matRad_calcSigmaIni(machine.data,stf(i).ray(j),stf(i).ray(j).SSD);
            
            if strcmp(pbCalcMode, 'fineSampling')
                % Ray tracing for beam i and ray j
                [ix,~,~,~,latDistsX,latDistsZ] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                                                     stf(i).sourcePoint_bev, ...
                                                     stf(i).ray(j).targetPoint_bev, ...
                                                     machine.meta.SAD, ...
                                                     find(~isnan(radDepthVdoseGrid{1})), ...
                                                     maxLateralCutoffDoseCalc);
                                                                  
                % Given the initial sigmas of the sampling ray, this
                % function provides the weights for the sub-pencil beams,
                % their positions and their sigma used for dose calculation
                for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
                    if (fineSamplingSigmaSub < sigmaIniRay(k)) && (fineSamplingSigmaSub > 0)
                        [finalWeight(:,k), sigmaSub(:,k), posX(:,k), posZ(:,k), numOfSub(:,k)] = ...
                                  matRad_calcWeights(sigmaIniRay(k), fineSamplingMethod, fineSamplingN, fineSamplingSigmaSub);
                    else
                        if (fineSamplingSigmaSub < 0)
                            matRad_cfg.dispError('Chosen fine sampling sigma cannot be negative!');
                        elseif (fineSamplingSigmaSub > sigmaIniRay(k))
                            matRad_cfg.dispError('Chosen fine sampling sigma is too high for defined plan!');
                        end                          
                    end
                end
            else
                % Ray tracing for beam i and ray j
                [ix,currRadialDist_sq,~,~,~,~] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                                                     stf(i).sourcePoint_bev, ...
                                                     stf(i).ray(j).targetPoint_bev, ...
                                                     machine.meta.SAD, ...
                                                     find(~isnan(radDepthVdoseGrid{1})), ...
                                                     maxLateralCutoffDoseCalc);
                                                                                  
                radDepths = radDepthVdoseGrid{1}(ix); 
            end
                   
            % just use tissue classes of voxels found by ray tracer
            if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                && strcmp(pln.radiationMode,'carbon')
                    vTissueIndex_j = vTissueIndex(ix,:);
            end
            
            

            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                counter = counter + 1;
                bixelsPerBeam = bixelsPerBeam + 1;
                
                % Display progress and update text only 200 times
                if mod(bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                        matRad_progress(bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                        floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                end
                
                % update waitbar only 100 times if it is not closed
                if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                    waitbar(counter/dij.totalNumOfBixels,figureWait);
                end

                % remember beam and bixel number
                if ~calcDoseDirect
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
                energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
                
                
                    if strcmp(pbCalcMode, 'fineSampling')
                    
                        % calculate projected coordinates for fine sampling of
                        % each beamlet
                        projCoords = matRad_projectOnComponents(VdoseGrid(ix), size(radDepthsMat{1}), stf(i).sourcePoint_bev,...
                                        stf(i).ray(j).targetPoint_bev, stf(i).isoCenter,...
                                        [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z],...
                                        -posX(:,k), -posZ(:,k), rotMat_system_T);

                        % interpolate radiological depths at projected
                        % coordinates
                        radDepths = interp3(radDepthsMat{1},projCoords(:,1,:)./dij.doseGrid.resolution.x,...
                            projCoords(:,2,:)./dij.doseGrid.resolution.y,projCoords(:,3,:)./dij.doseGrid.resolution.z,'nearest');                       

                        % compute radial distances relative to pencil beam
                        % component
                        currRadialDist_sq = reshape(bsxfun(@plus,latDistsX,posX(:,k)'),[],1,numOfSub(k)).^2 + reshape(bsxfun(@plus,latDistsZ,posZ(:,k)'),[],1,numOfSub(k)).^2;
                    end
                
                    % create offset vector to account for additional offsets modelled in the base data and a potential 
                    % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                    % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow  
                    % interpolations in matRad_calcParticleDoseBixel() and avoid extrapolations.
                    offsetRadDepth = machine.data(energyIx).offset - stf(i).ray(j).rangeShifter(k).eqThickness;

                    % find depth depended lateral cut off
                    if cutOffLevel >= 1
                        currIx = radDepths <= machine.data(energyIx).depths(end) + offsetRadDepth;
                    elseif cutOffLevel < 1 && cutOffLevel > 0
                        % perform rough 2D clipping
                        currIx = radDepths <= machine.data(energyIx).depths(end) + offsetRadDepth & ...
                             currRadialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);

                        % peform fine 2D clipping  
                        if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                            currIx(currIx) = matRad_interp1((machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                (machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= currRadialDist_sq(currIx);
                        end
                    else
                        matRad_cfg.dispError('Cutoff must be a value between 0 and 1!')
                    end

                    % empty bixels may happen during recalculation of error
                    % scenarios -> skip to next bixel
                    if ~(any(any(currIx)))
                        continue;
                    end

                    % adjust radDepth according to range shifter
                    currRadDepths = radDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness;

                    % select correct initial focus sigma squared
                    sigmaIni_sq = sigmaIniRay(k)^2;
                          
                    % consider range shifter for protons if applicable
                    if stf(i).ray(j).rangeShifter(k).eqThickness > 0 && strcmp(pln.radiationMode,'protons')

                        % compute!
                        sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx).energy, ...
                                                           stf(i).ray(j).rangeShifter(k), ...
                                                           stf(i).ray(j).SSD);

                        % add to initial sigma in quadrature
                        sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;
                        
                    end
                                  
                if strcmp(pbCalcMode, 'fineSampling')
                    % initialise empty dose array
                    totalDose = zeros(size(currIx,1),1);
                    
                    if isfield(dij,'mLETDose')
                        % calculate particle LET for bixel k on ray j of beam i
                        depths = machine.data(energyIx).depths + machine.data(energyIx).offset; 
                        totalLET = zeros(size(currIx,1),1);
                    end
                    
                    % run over components
                    for c = 1:numOfSub(k)
                        tmpDose = zeros(size(currIx,1),1);
                        bixelDose = finalWeight(c,k).*matRad_calcParticleDoseBixel(...
                                radDepths(currIx(:,:,c),1,c), ...
                                currRadialDist_sq(currIx(:,:,c),:,c), ...
                                sigmaSub(k)^2, ...
                                machine.data(energyIx));
                                                        
                        tmpDose(currIx(:,:,c)) = bixelDose;
                        totalDose = totalDose + tmpDose;
                        
                        if isfield(dij,'mLETDose') 
                            tmpLET = zeros(size(currIx,1),1);
                            tmpLET(currIx(:,:,c)) = matRad_interp1(depths,machine.data(energyIx).LET,radDepths(currIx(:,:,c),1,c));    
                            totalLET = totalLET + tmpLET;
                        end
                    end
                    
                    doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(VdoseGrid(ix),1,totalDose,dij.doseGrid.numOfVoxels,1);
                    if isfield(dij,'mLETDose') 
                        letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(VdoseGrid(ix),1,totalDose.*totalLET,dij.doseGrid.numOfVoxels,1);
                    end                    
                else
                    % calculate particle dose for bixel k on ray j of beam i
                    bixelDose = matRad_calcParticleDoseBixel(...
                        currRadDepths, ...
                        currRadialDist_sq(currIx), ...
                        sigmaIni_sq, ...
                        machine.data(energyIx));                 

                    % dij sampling is exluded for particles until we investigated the influence of voxel sampling for particles
                    %relDoseThreshold   =  0.02;   % sample dose values beyond the relative dose
                    %Type               = 'dose';
                    %[currIx,bixelDose] = matRad_DijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);

                    % Save dose for every bixel in cell array
                    doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(VdoseGrid(ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);
                
                    if isfield(dij,'mLETDose')
                      % calculate particle LET for bixel k on ray j of beam i
                      depths = machine.data(energyIx).depths + machine.data(energyIx).offset; 
                      bixelLET = matRad_interp1(depths,machine.data(energyIx).LET,currRadDepths); 

                      % Save LET for every bixel in cell array
                      letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(VdoseGrid(ix(currIx)),1,bixelLET.*bixelDose,dij.doseGrid.numOfVoxels,1);
                    end
                end

                
                             
                if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                    && strcmp(pln.radiationMode,'carbon')
                    % calculate alpha and beta values for bixel k on ray j of                  
                    [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                        currRadDepths,...
                        vTissueIndex_j(currIx,:),...
                        machine.data(energyIx));
                    
                    alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(VdoseGrid(ix(currIx)),1,bixelAlpha.*bixelDose,dij.doseGrid.numOfVoxels,1);
                    betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1}  = sparse(VdoseGrid(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.doseGrid.numOfVoxels,1);
                end
                
                matRad_calcDoseFillDij;                

            end
            
        end
        
    end
end

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end

    
