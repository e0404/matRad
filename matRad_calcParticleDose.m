function dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle dose calculation wrapper
% 
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
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

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln);

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for particles...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if calcDoseDirect 
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
end

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
doseTmpContainer     = cell(numOfBixelsContainer,dij.numOfScenarios);

if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
        && strcmp(pln.radiationMode,'carbon')
   
        alphaDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
        betaDoseTmpContainer  = cell(numOfBixelsContainer,dij.numOfScenarios);
        for i = 1:dij.numOfScenarios
            dij.mAlphaDose{i}    = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
            dij.mSqrtBetaDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
        end
        
elseif isequal(pln.propOpt.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
            dij.RBE = 1.1;
            fprintf(['matRad: Using a constant RBE of 1.1 \n']);   
end

% Only take voxels inside patient.
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(dij.numOfVoxels,1);
eraseCtDensMask(V) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

if isfield(pln,'propDoseCalc') && ...
   isfield(pln.propDoseCalc,'calcLET') && ...
   pln.propDoseCalc.calcLET
  if isfield(machine.data,'LET')
    letDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
    % Allocate space for dij.dosexLET sparse matrix
    for i = 1:dij.numOfScenarios
        dij.mLETDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
    end
  else
    warndlg('LET not available in the machine data. LET will not be calculated.');
  end
end

% generates tissue class matrix for biological optimization
if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
        && strcmp(pln.radiationMode,'carbon')
    
    fprintf('matRad: loading biological base data... ');
    vTissueIndex = zeros(size(V,1),1);
    
    %set overlap priorities
    cst  = matRad_setOverlapPriorities(cst);
    
    for i = 1:size(cst,1)
        % find indices of structures related to V
        [~, row] = ismember(vertcat(cst{i,4}{:}),V,'rows'); 
        % check if base data contains alphaX and betaX
        if   isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
            % check if cst is compatiable 
            if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX') 

                IdxTissue = find(ismember(machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                                 ismember(machine.data(1).betaX,cst{i,5}.betaX));

                % check consitency of biological baseData and cst settings
                if ~isempty(IdxTissue)
                    vTissueIndex(row) = IdxTissue;
                else
                    error('biological base data and cst inconsistent\n');
                end
            else
                vTissueIndex(row) = 1;
                fprintf(['matRad: tissue type of ' cst{i,2} ' was set to 1 \n']);
            end
        else
            error('base data is incomplement - alphaX and/or betaX is missing');
        end
        
    end
    fprintf('done.\n');

% issue warning if biological optimization not possible
elseif sum(strcmp(pln.propOpt.bioOptimization,{'LEMIV_effect','LEMIV_RBExD'}))>0 && ~strcmp(pln.radiationMode,'carbon') ||...
       ~strcmp(pln.radiationMode,'protons') && strcmp(pln.propOpt.bioOptimization,'const_RBExD')
    warndlg([pln.propOpt.bioOptimization ' optimization not possible with ' pln.radiationMode '- physical optimization is carried out instead.']);
    pln.propOpt.bioOptimization = 'none';      
end

% compute SSDs
stf = matRad_computeSSD(stf,ct);

fprintf('matRad: Particle dose calculation...\n');
counter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(stf) % loop over all beams
    
    fprintf(['Beam ' num2str(i) ' of ' num2str(length(stf)) ': \n']);
        
    % remember beam and bixel number
    if calcDoseDirect
        dij.beamNum(i)    = i;
        dij.rayNum(i)     = i;
        dij.bixelNum(i)   = i;
    end

    bixelsPerBeam = 0;
    
    % convert voxel indices to real coordinates using iso center of beam i
    xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
    yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
    zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
    coordsV  = [xCoordsV yCoordsV zCoordsV];

    % Get Rotation Matrix
    % Do not transpose matrix since we usage of row vectors &
    % transformation of the coordinate system need double transpose

    % rotation around Z axis (gantry)
    rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

    % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
    rot_coordsV = coordsV*rotMat_system_T;

    rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

    % Calcualte radiological depth cube
    lateralCutoffRayTracing = 50;
    fprintf('matRad: calculate radiological depth cube...');
    radDepthV = matRad_rayTracing(stf(i),ct,V,rot_coordsV,lateralCutoffRayTracing);
    fprintf('done.\n');
    
    % get indices of voxels where ray tracing results are available
    radDepthIx = find(~isnan(radDepthV{1}));
    
    % limit rotated coordinates to positions where ray tracing is availabe
    rot_coordsV = rot_coordsV(radDepthIx,:);
    
    % Determine lateral cutoff
    fprintf('matRad: calculate lateral cutoff...');
    cutOffLevel = .995;
    visBoolLateralCutOff = 0;
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
    fprintf('done.\n');    

    for j = 1:stf(i).numOfRays % loop over all rays

        if ~isempty(stf(i).ray(j).energy)

            % find index of maximum used energy (round to keV for numerical
            % reasons
            energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);

            maxLateralCutoffDoseCalc = max(machine.data(energyIx).LatCutOff.CutOff);

            % Ray tracing for beam i and ray j
            [ix,radialDist_sq] = matRad_calcGeoDists(rot_coordsV, ...
                                                     stf(i).sourcePoint_bev, ...
                                                     stf(i).ray(j).targetPoint_bev, ...
                                                     machine.meta.SAD, ...
                                                     radDepthIx, ...
                                                     maxLateralCutoffDoseCalc);
            radDepths = radDepthV{1}(ix);   
            
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
                end
                
                % find energy index in base data
                energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
                
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
                         radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);

                    % peform fine 2D clipping  
                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                        currIx(currIx) = matRad_interp1((machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                            (machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= radialDist_sq(currIx);
                    end
                else
                    error('cutoff must be a value between 0 and 1')
                end
                
                % empty bixels may happen during recalculation of error
                % scenarios -> skip to next bixel
                if ~any(currIx)
                    continue;
                end
                
                % adjust radDepth according to range shifter
                currRadDepths = radDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness;

                % calculate initial focus sigma
                sigmaIni = matRad_interp1(machine.data(energyIx).initFocus.dist (stf(i).ray(j).focusIx(k),:)', ...
                                             machine.data(energyIx).initFocus.sigma(stf(i).ray(j).focusIx(k),:)',stf(i).ray(j).SSD);
                sigmaIni_sq = sigmaIni^2;
                
                % consider range shifter for protons if applicable
                if stf(i).ray(j).rangeShifter(k).eqThickness > 0 && strcmp(pln.radiationMode,'protons')
                    
                    % compute!
                    sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx).energy, ...
                                                       stf(i).ray(j).rangeShifter(k), ...
                                                       stf(i).ray(j).SSD);
                              
                    % add to initial sigma in quadrature
                    sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;
                    
                end
                                
                % calculate particle dose for bixel k on ray j of beam i
                bixelDose = matRad_calcParticleDoseBixel(...
                    currRadDepths, ...
                    radialDist_sq(currIx), ...
                    sigmaIni_sq, ...
                    machine.data(energyIx));                 
  
                % dij sampling is exluded for particles until we investigated the influence of voxel sampling for particles
                %relDoseThreshold   =  0.02;   % sample dose values beyond the relative dose
                %Type               = 'dose';
                %[currIx,bixelDose] = matRad_DijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);
                
                % Save dose for every bixel in cell array
                doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelDose,dij.numOfVoxels,1);

                if isfield(dij,'mLETDose')
                  % calculate particle LET for bixel k on ray j of beam i
                  depths = machine.data(energyIx).depths + machine.data(energyIx).offset; 
                  bixelLET = matRad_interp1(depths,machine.data(energyIx).LET,radDepths(currIx)); 

                  % Save LET for every bixel in cell array
                  letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelLET.*bixelDose,dij.numOfVoxels,1);
                end
                             
                if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                    && strcmp(pln.radiationMode,'carbon')
                    % calculate alpha and beta values for bixel k on ray j of                  
                    [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                        radDepths(currIx),...
                        vTissueIndex_j(currIx,:),...
                        machine.data(energyIx));
                
                    alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelAlpha.*bixelDose,dij.numOfVoxels,1);
                    betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1}  = sparse(V(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.numOfVoxels,1);
                end
                
                % save computation time and memory by sequentially filling the
                % sparse matrix dose.dij from the cell array
                if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                    
                    if calcDoseDirect
                        if isfield(stf(1).ray(1),'weight') && numel(stf(i).ray(j).weight) >= k
                            
                            % score physical dose
                            dij.physicalDose{1}(:,i) = dij.physicalDose{1}(:,i) + stf(i).ray(j).weight(k) * doseTmpContainer{1,1};
                            
                            if isfield(dij,'mLETDose')
                                dij.mLETDose{1}(:,i) = dij.mLETDose{1}(:,i) + stf(i).ray(j).weight(k) * letDoseTmpContainer{1,1};
                            end
                            
                            if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                                && strcmp(pln.radiationMode,'carbon')
                                
                                % score alpha and beta matrices
                                dij.mAlphaDose{1}(:,i)    = dij.mAlphaDose{1}(:,i) + stf(i).ray(j).weight(k) * alphaDoseTmpContainer{1,1};
                                dij.mSqrtBetaDose{1}(:,i) = dij.mSqrtBetaDose{1}(:,i) + stf(i).ray(j).weight(k) * betaDoseTmpContainer{1,1};
                                
                            end
                        else
                            
                            error(['No weight available for beam ' num2str(i) ', ray ' num2str(j) ', bixel ' num2str(k)]);
                            
                        end
                    else
                        
                        dij.physicalDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                        
                        if isfield(dij,'mLETDose')
                            dij.mLETDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [letDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                        end
                        
                        if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                            && strcmp(pln.radiationMode,'carbon')
                        
                            dij.mAlphaDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [alphaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                            dij.mSqrtBetaDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [betaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                        end
                    
                    end
                end
                

            end
            
        end
        
    end
end

try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end

    