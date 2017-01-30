function dij = matRad_calcParticleDose(ct,stf,pln,cst,multScen,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle dose calculation wrapper
% cst
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   multScen:       matRad multiple scnerio struct
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

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for particles...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');

% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;
dij.numOfScenarios     = multScen.totalNumOfScen;
dij.ScenProb           = multScen.shiftScenProb;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
for CtScen = 1:multScen.numOfCtScen
    for ShiftScen = 1:multScen.numOfShiftScen
        for RangeShiftScen = 1:multScen.numOfRangeShiftScen  
            
            if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)
                dij.physicalDose{CtScen,ShiftScen,RangeShiftScen} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
            end
            
        end
    end
end

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
if calcDoseDirect
    numOfBixelsContainer = 1;
else
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

doseTmpContainer = cell(numOfBixelsContainer,multScen.numOfCtScen,multScen.numOfShiftScen,multScen.numOfRangeShiftScen);
if (isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) && strcmp(pln.radiationMode,'carbon') || ...
   (isequal(pln.bioOptimization,'LSM_effect')   || isequal(pln.bioOptimization,'LSM_RBExD'))   && strcmp(pln.radiationMode,'protons')

    alphaDoseTmpContainer = cell(numOfBixelsContainer,multScen.numOfCtScen,multScen.numOfShiftScen,multScen.numOfRangeShiftScen);
    betaDoseTmpContainer  = cell(numOfBixelsContainer,multScen.numOfCtScen,multScen.numOfShiftScen,multScen.numOfRangeShiftScen);
    
    for CtScen = 1:multScen.numOfCtScen
        for ShiftScen = 1:multScen.numOfShiftScen
            for RangeShiftScen = 1:multScen.numOfRangeShiftScen  
            
                if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)
                    dij.mAlphaDose{CtScen,ShiftScen,RangeShiftScen}        = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
                    dij.mSqrtBetaDose{CtScen,ShiftScen,RangeShiftScen}     = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
                end
                
            end
        end
    end
    
    % define additional parameter for the LSM model
    if (isequal(pln.bioOptimization,'LSM_effect')   || isequal(pln.bioOptimization,'LSM_RBExD'))   && strcmp(pln.radiationMode,'protons')
           pln.calcLET            = true;
           dij.lamda_1_1          = 0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (fitted for head and neck patients)
           dij.corrFacEntranceRBE = 0.5;   %[kev/mum]
           dij.upperLETThreshold  = 30;    %[kev/mum]
           dij.lowerLETThreshold  = 0.3;    %[kev/mum]
    end

elseif isequal(pln.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
    dij.RBE = 1.1;
    fprintf(['matRad: Using a constant RBE of 1.1 \n']);   
end

% Only take voxels inside patient.
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

 % allocate space for dij.dosexLET sparse matrix
if (isfield(pln,'calcLET') && pln.calcLET) || (isequal(pln.bioOptimization,'LSM_effect') || isequal(pln.bioOptimization,'LSM_RBExD')) && strcmp(pln.radiationMode,'protons')
  if isfield(machine.data,'LET')
      
    letDoseTmpContainer = cell(numOfBixelsContainer,multScen.numOfCtScen,multScen.numOfShiftScen,multScen.numOfRangeShiftScen);
   
    for CtScen = 1:multScen.numOfCtScen
        for ShiftScen = 1:multScen.numOfShiftScen
            for RangeShiftScen = 1:multScen.numOfRangeShiftScen  
            
                if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)
                     dij.mLETDose{CtScen,ShiftScen,RangeShiftScen} = spalloc(prod(ct.cubeDim),dij.totalNumOfBixels,1);
                end
                
            end
        end
    end
  else
    warndlg('LET not available in the machine data. LET will not be calculated.');
  end
end


    
% generate two vectors depicting reference radiosensitivity parameter of photons
if (isequal(pln.bioOptimization,'LSM_effect')   || isequal(pln.bioOptimization,'LSM_RBExD')) && strcmp(pln.radiationMode,'protons')

    fprintf('matRad: precumputations for biological treatment planning... ');
    alphaX = zeros(size(V,1),1);
    betaX  = zeros(size(V,1),1);
    %set overlap priorities
    cst  = matRad_setOverlapPriorities(cst);
    
    for i = 1:size(cst,1)
        % find indices of structures related to V
        [~, row] = ismember(vertcat(cst{i,4}{:}),V,'rows'); 
       
        % check if cst is compatiable 
        if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX') 
            alphaX(row) = cst{i,5}.alphaX;
            betaX(row)  = cst{i,5}.betaX;               
        else
            alphaX(row) = 0.1;    % default parameter  
            betaX(row)  = 0.05;   % default parameter
            fprintf(['matRad: using default alpha_x and beta_x parameters for ' cst{i,2} ' \n']);
        end
  
    end
    fprintf('done.\n');
        
% generates tissue class matrix for biological optimization
elseif (isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) ... 
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
elseif sum(strcmp(pln.bioOptimization,{'LEMIV_effect','LEMIV_RBExD'}))>0 && ~strcmp(pln.radiationMode,'carbon') ||...
       ~strcmp(pln.radiationMode,'protons') && strcmp(pln.bioOptimization,'const_RBExD')
    warndlg([pln.bioOptimization ' optimization not possible with ' pln.radiationMode '- physical optimization is carried out instead.']);
    pln.bioOptimization = 'none';      
end

for ShiftScen = 1:multScen.numOfShiftScen

    % manipulate isocenter
    pln.isoCenter    = pln.isoCenter + multScen.shifts(:,ShiftScen)';
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + multScen.shifts(:,ShiftScen)';
    end
    
    fprintf(['shift scenario ' num2str(ShiftScen) ' of ' num2str(multScen.numOfShiftScen) ': \n']);
    fprintf('matRad: Particle dose calculation...\n');
    counter = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:dij.numOfBeams % loop over all beams

        fprintf(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ': \n']);

        bixelsPerBeam = 0;

        % convert voxel indices to real coordinates using iso center of beam i
        xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
        yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
        zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
        coordsV  = [xCoordsV yCoordsV zCoordsV];

        % Set gantry and couch rotation matrices according to IEC 61217
        % Use transpose matrices because we are working with row vectros

        % rotation around Z axis (gantry)
        inv_rotMx_XY_T = [ cosd(-pln.gantryAngles(i)) sind(-pln.gantryAngles(i)) 0;
                          -sind(-pln.gantryAngles(i)) cosd(-pln.gantryAngles(i)) 0;
                                                    0                          0 1];
        % rotation around Y axis (couch)
        inv_rotMx_XZ_T = [cosd(-pln.couchAngles(i)) 0 -sind(-pln.couchAngles(i));
                                                  0 1                         0;
                          sind(-pln.couchAngles(i)) 0  cosd(-pln.couchAngles(i))];

        % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
        rot_coordsV = coordsV*inv_rotMx_XZ_T*inv_rotMx_XY_T;
        
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
        cutOffLevel = .99;
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

                % just use tissue classes of voxels found by ray tracer
                if (isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) ... 
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

                    % remember beam and  bixel number
                    dij.beamNum(counter)  = i;
                    dij.rayNum(counter)   = j;
                    dij.bixelNum(counter) = k;

                    % find energy index in base data
                    energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));

                    for CtScen = 1:multScen.numOfCtScen
                        for RangeShiftScen = 1:multScen.numOfRangeShiftScen  

                            if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)

                                % manipulate radDepthCube for range scenarios
                                radDepths = radDepthV{CtScen}(ix);                                          

                                if multScen.relRangeShifts(RangeShiftScen) ~= 0 || multScen.absRangeShifts(RangeShiftScen) ~= 0
                                    radDepths = radDepths +...                                                                                % original cube
                                                radDepthV{CtScen}(ix)*multScen.relRangeShifts(RangeShiftScen) +... % rel range shift
                                                multScen.absRangeShifts(RangeShiftScen);                                                      % absolute range shift
                                    radDepths(radDepths < 0) = 0;  
                                end

                                % find depth depended lateral cut off
                                if cutOffLevel >= 1
                                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset;
                                elseif cutOffLevel < 1 && cutOffLevel > 0
                                    % perform rough 2D clipping
                                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset & ...
                                         radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);

                                    % peform fine 2D clipping  
                                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                                        currIx(currIx) = matRad_interp1((machine.data(energyIx).LatCutOff.depths + machine.data(energyIx).offset)',...
                                            (machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= radialDist_sq(currIx);
                                    end
                                else
                                    error('cutoff must be a value between 0 and 1')
                                end

                                % calculate particle dose for bixel k on ray j of beam i
                                bixelDose = matRad_calcParticleDoseBixel(...
                                    radDepths(currIx), ...
                                    radialDist_sq(currIx), ...
                                    stf(i).ray(j).SSD{CtScen}, ...
                                    stf(i).ray(j).focusIx(k), ...
                                    machine.data(energyIx)); 


                                % dij sampling is exluded for particles until we investigated the influence of voxel sampling for particles
                                %relDoseThreshold   =  0.02;   % sample dose values beyond the relative dose
                                %Type               = 'dose';
                                %[currIx,bixelDose] = matRad_DijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);


                                % Save dose for every bixel in cell array
                                doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen} = sparse(V(ix(currIx)),1,bixelDose,dij.numOfVoxels,1); 

                                if isfield(dij,'mLETDose')
                                  % calculate particle LET for bixel k on ray j of beam i
                                  depths = machine.data(energyIx).depths + machine.data(energyIx).offset; 
                                  bixelLET = matRad_interp1(depths,machine.data(energyIx).LET,radDepths(currIx)); 

                                  % Save LET for every bixel in cell array
                                  letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen} = sparse(V(ix(currIx)),1,bixelLET.*bixelDose,dij.numOfVoxels,1);
                                end                            

                                if ((isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) && strcmp(pln.radiationMode,'carbon')) ||...
                                   ((isequal(pln.bioOptimization,'LSM_effect')   || isequal(pln.bioOptimization,'LSM_RBExD'))   && strcmp(pln.radiationMode,'protons'))

                                     if ((isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD'))  && strcmp(pln.radiationMode,'carbon'))

                                         % calculate alpha and beta values for bixel k on ray j of                  
                                        [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                                            radDepths(currIx),...
                                            vTissueIndex_j(currIx,:),...
                                            machine.data(energyIx));

                                     elseif ((isequal(pln.bioOptimization,'LSM_effect')   || isequal(pln.bioOptimization,'LSM_RBExD'))   && strcmp(pln.radiationMode,'protons'))

                                        alpha_0            = alphaX(ix(currIx)) - (dij.lamda_1_1 * dij.corrFacEntranceRBE);  
                                        bixelAlpha         = zeros(size(bixelDose)); bixelBeta = zeros(size(bixelDose));
                                        ixLSM              = dij.lowerLETThreshold < bixelLET < dij.upperLETThreshold;

                                        bixelAlpha(ixLSM)  = alpha_0(ixLSM) + dij.lamda_1_1 * bixelLET;
                                        
                                        if sum(ixLSM) < length(bixelLET)
                                            bixelAlpha(~ixLSM)                           = (alpha_0(~ixLSM) + dij.lamda_1_1 * dij.upperLetThreshold)./(dij.upperLetThreshold/bixelLET);
                                            bixelAlpha(bixelLET < dij.lowerLETThreshold) =  alpha_0(bixelLET < dij.lowerLETThreshold) + dij.lamda_1_1 * dij.lowerLETThreshold;
                                        end
                                        
                                        bixelBeta(:)         = betaX(ix(currIx));

                                     end

                                    alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen} = sparse(V(ix(currIx)),1,bixelAlpha.*bixelDose,dij.numOfVoxels,1);
                                    betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen}  = sparse(V(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.numOfVoxels,1);
                                end
                                
                            end

                        end
                    end

                    % save computation time and memory by sequentially filling the
                    % sparse matrix dose.dij from the cell array
                    if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                        for CtScen = 1:multScen.numOfCtScen
                            for RangeShiftScen = 1:multScen.numOfRangeShiftScen

                                if multScen.ScenCombMask(CtScen,ShiftScen,RangeShiftScen)
                                    if calcDoseDirect
                                        if isfield(stf(1).ray(1),'weight') && numel(stf(i).ray(j).weight) >= k    
                                            % score physical dose
                                            dij.physicalDose{CtScen,ShiftScen,RangeShiftScen}(:,1) = dij.physicalDose{CtScen,ShiftScen,RangeShiftScen}(:,1) + stf(i).ray(j).weight(k) * doseTmpContainer{1,CtScen,ShiftScen,RangeShiftScen};

                                            if isfield(dij,'mLETDose')
                                                % score LETxDose matrices
                                                dij.mLETDose{CtScen,ShiftScen,RangeShiftScen}(:,1) = dij.mLETDose{CtScen,ShiftScen,RangeShiftScen}(:,1) + stf(i).ray(j).weight(k) * letDoseTmpContainer{1,CtScen,ShiftScen,RangeShiftScen}; 
                                            end
                                            if (isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) && strcmp(pln.radiationMode,'carbon') || ...
                                               ((isequal(pln.bioOptimization,'LSM_effect')  || isequal(pln.bioOptimization,'LSM_RBExD'))   && strcmp(pln.radiationMode,'protons'))

                                                % score alphaxDose and sqrt(beta)xDose matrices
                                                dij.mAlphaDose{CtScen,ShiftScen,RangeShiftScen}(:,1)    = dij.mAlphaDose{CtScen,ShiftScen,RangeShiftScen}(:,1)    + stf(i).ray(j).weight(k) * alphaDoseTmpContainer{1,CtScen,ShiftScen,RangeShiftScen};
                                                dij.mSqrtBetaDose{CtScen,ShiftScen,RangeShiftScen}(:,1) = dij.mSqrtBetaDose{CtScen,ShiftScen,RangeShiftScen}(:,1) + stf(i).ray(j).weight(k) * betaDoseTmpContainer{1,CtScen,ShiftScen,RangeShiftScen};
                                            end
                            
                                        else
                                            error(['No weight available for beam ' num2str(i) ', ray ' num2str(j) ', bixel ' num2str(k)]);
                                        end
                                    else
                                        % fill entire dose influence matrix
                                        dij.physicalDose{CtScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen}];

                                        if isfield(dij,'mLETDose')
                                            % fill entire LETxDose influence matrix
                                            dij.mLETDose{CtScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [letDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen}]; 
                                        end
                                        
                                        if ((isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) && strcmp(pln.radiationMode,'carbon')) || ...
                                           ((isequal(pln.bioOptimization,'LSM_effect')   || isequal(pln.bioOptimization,'LSM_RBExD'))   && strcmp(pln.radiationMode,'protons'))
                                            % fill entire alphaxDose influence and sqrt(beta)xDose influence matrices
                                            dij.mAlphaDose{CtScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter)    = [alphaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen}];
                                            dij.mSqrtBetaDose{CtScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [betaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,CtScen,ShiftScen,RangeShiftScen}];
                                        end
                                        
                                    end

                                end

                            end
                            
                        end

                    end

                end

            end

        end
        
    end

    % manipulate isocenter
    pln.isoCenter    = pln.isoCenter - multScen.shifts(:,ShiftScen)';
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter - multScen.shifts(:,ShiftScen)';
    end 

end

% set dummy value
dij.indexforOpt = 1;

try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end
