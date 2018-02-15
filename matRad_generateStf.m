function stf = matRad_generateStf(ct,cst,pln,visMode)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad steering information generation
%
% call
%   stf = matRad_generateStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
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


fprintf('matRad: Generating stf struct... ');

if nargin < 4
    visMode = 0;
end

if numel(pln.propStf.gantryAngles) ~= numel(pln.propStf.couchAngles)
    error('Inconsistent number of gantry and couch angles.');
end

if pln.propStf.bixelWidth < 0 || ~isfinite(pln.propStf.bixelWidth)
   error('bixel width (spot distance) needs to be a real number [mm] larger than zero.');
end

% find all target voxels from cst cell array
V = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;vertcat(cst{i,4}{:})];
    end
end

% Remove double voxels
V = unique(V);

% generate voi cube for targets
voiTarget    = zeros(ct.cubeDim);
voiTarget(V) = 1;

% add margin
addmarginBool = 0;
if addmarginBool
    voiTarget = matRad_addMargin(voiTarget,cst,ct.resolution,ct.resolution,true);
    V   = find(voiTarget>0);
end

% throw error message if no target is found
if isempty(V)
    error('Could not find target.');
end

% prepare structures necessary for particles
fileName = [pln.radiationMode '_' pln.machine];
try
    load([fileparts(mfilename('fullpath')) filesep fileName]);
    SAD = machine.meta.SAD;
catch
    error(['Could not find the following machine file: ' fileName ]);
end

if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    
    availableEnergies = [machine.data.energy];
    availablePeakPos  = [machine.data.peakPos] + [machine.data.offset];
    
    if sum(availablePeakPos<0)>0
        error('at least one available peak position is negative - inconsistent machine file')
    end
    %clear machine;
end

% Convert linear indices to 3D voxel coordinates
[coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(ct.cubeDim,V);

% Define steering file like struct. Prellocating for speed.
% Preallocate with known size of array to avoid errors in VMAT-specific
% code
stf = struct('gantryAngle',cell(size(pln.propStf.gantryAngles)));

if pln.propOpt.runVMAT
    %Initialize master ray positions and target points with NaNs, to be
    %deleted later.  These arrays are the unions of the corresponding
    %arrays per gantry angle.  In order to do VMAT, it is easier to have
    %the same MLC range and dij calculation for every possible beam/gantry
    %angle.
    masterRayPosBEV = nan(1,3);
    masterTargetPointBEV = nan(1,3);
end


% loop over all angles
for i = 1:length(pln.propStf.gantryAngles)
    
    % Correct for iso center position. Whit this correction Isocenter is
    % (0,0,0) [mm]
    coordsX = coordsX_vox*ct.resolution.x - pln.propStf.isoCenter(i,1);
    coordsY = coordsY_vox*ct.resolution.y - pln.propStf.isoCenter(i,2);
    coordsZ = coordsZ_vox*ct.resolution.z - pln.propStf.isoCenter(i,3);
    
    % Save meta information for treatment plan
    stf(i).gantryAngle   = pln.propStf.gantryAngles(i);
    stf(i).couchAngle    = pln.propStf.couchAngles(i);
    stf(i).bixelWidth    = pln.propStf.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    stf(i).SAD           = SAD;
    stf(i).isoCenter     = pln.propStf.isoCenter(i,:);
    
    % Get the (active) rotation matrix. We perform a passive/system
    % rotation with row vector coordinates, which would introduce two
    % inversions / transpositions of the matrix, thus no changes to the
    % rotation matrix are necessary
    rotMat_system_T = matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i));
    
    rot_coords = [coordsX coordsY coordsZ]*rotMat_system_T;
    
    % project x and z coordinates to isocenter
    coordsAtIsoCenterPlane(:,1) = (rot_coords(:,1)*SAD)./(SAD + rot_coords(:,2));
    coordsAtIsoCenterPlane(:,2) = (rot_coords(:,3)*SAD)./(SAD + rot_coords(:,2));
    
    % Take unique rows values for beamlets positions. Calculate position of
    % central ray for every bixel    
    rayPos = unique(pln.propStf.bixelWidth*round([           coordsAtIsoCenterPlane(:,1) ... 
                                                  zeros(size(coordsAtIsoCenterPlane,1),1) ...
                                                             coordsAtIsoCenterPlane(:,2)]/pln.propStf.bixelWidth),'rows');
                                                  
    % pad ray position array if resolution of target voxel grid not sufficient
    maxCtResolution = max([ct.resolution.x ct.resolution.y ct.resolution.z]);
    if pln.propStf.bixelWidth < maxCtResolution
        origRayPos = rayPos;
        for j = -floor(maxCtResolution/pln.propStf.bixelWidth):floor(maxCtResolution/pln.propStf.bixelWidth)
            for k = -floor(maxCtResolution/pln.propStf.bixelWidth):floor(maxCtResolution/plnpropStf.bixelWidth)
                if abs(j)+abs(k)==0
                    continue;
                end                
                rayPos = [rayPos; origRayPos(:,1)+j*pln.propStf.bixelWidth origRayPos(:,2) origRayPos(:,3)+k*pln.propStf.bixelWidth];
            end
        end
     end

     % remove spaces within rows of bixels for DAO
     if pln.propOpt.runDAO
         % create single x,y,z vectors
         x = rayPos(:,1);
         y = rayPos(:,2);
         z = rayPos(:,3);
         uniZ = unique(z);
         for j = 1:numel(uniZ)
             x_loc = x(z == uniZ(j));
             x_min = min(x_loc);
             x_max = max(x_loc);
             x = [x; [x_min:pln.propStf.bixelWidth:x_max]'];
             y = [y; zeros((x_max-x_min)/pln.propStf.bixelWidth+1,1)];
             z = [z; uniZ(j)*ones((x_max-x_min)/pln.propStf.bixelWidth+1,1)];             
         end
         
         rayPos = [x,y,z];
     end
    
    % remove double rays
    rayPos = unique(rayPos,'rows');
    
    % Save the number of rays
    stf(i).numOfRays = size(rayPos,1);
    
    % Save ray and target position in beam eye's view (bev)
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos_bev = rayPos(j,:);
        stf(i).ray(j).targetPoint_bev = [2*stf(i).ray(j).rayPos_bev(1) ...
            SAD ...
            2*stf(i).ray(j).rayPos_bev(3)];
    end
    
    if ~pln.propOpt.runVMAT
        %If it is VMAT, we will do this later
        
        % source position in bev
        stf(i).sourcePoint_bev = [0 -SAD 0];
        
        % get (active) rotation matrix
        % transpose matrix because we are working with row vectors
        rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i)));
        
        
        stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;
        
        % Save ray and target position in lps system.
        for j = 1:stf(i).numOfRays
            stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMat_vectors_T;
            stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMat_vectors_T;
            if strcmp(pln.radiationMode,'photons')
                stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                    [rayPos(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    rayPos(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    rayPos(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                    rayPos(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMat_vectors_T;
            end
        end
        
        % loop over all rays to determine meta information for each ray
        stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
        
        
        for j = stf(i).numOfRays:-1:1
            
            % ray tracing necessary to determine depth of the target
            [~,l,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                ct.resolution, ...
                stf(i).sourcePoint, ...
                stf(i).ray(j).targetPoint, ...
                [ct.cube {voiTarget}]);
                
                % find appropriate energies for particles
                if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')
                    
                    % target hit
                    if sum(rho{2}) > 0
                        
                        % compute radiological depths
                        % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                        radDepths = cumsum(l .* rho{1});
                        
                        % find target entry & exit
                        diff_voi    = diff([rho{2}]);
                        targetEntry = radDepths(diff_voi == 1);
                        targetExit  = radDepths(diff_voi == -1);
                        
                        if numel(targetEntry) ~= numel(targetExit)
                            error('Inconsistency during ray tracing.');
                        end
                        
                        stf(i).ray(j).energy = [];
                        
                        % Save energies in stf struct
                        for k = 1:numel(targetEntry)
                            stf(i).ray(j).energy = [stf(i).ray(j).energy availableEnergies(availablePeakPos>=targetEntry(k)&availablePeakPos<=targetExit(k))];
                            % adjust spot spacing according to pln.bixelWidth when using HIT basedata
                            %DefaultLongitudialSpotSpacing = pln.bixelWidth;  % in [mm]
                            DefaultLongitudialSpotSpacing = 3;
                            if strcmp(pln.machine,'HIT') && length(stf(i).ray(j).energy)>2
                                Tolerance = 0.5;
                                hasNext = true;
                                CntEnergy =2;
                                while hasNext
                                    if abs(stf(i).ray(j).energy(CntEnergy)-stf(i).ray(j).energy(CntEnergy-1))<...
                                            DefaultLongitudialSpotSpacing-Tolerance
                                        stf(i).ray(j).energy(CntEnergy)=[];
                                    else
                                        CntEnergy = CntEnergy+1;
                                    end
                                    if CntEnergy == length(stf(i).ray(j).energy)
                                        hasNext = false;
                                    end
                                end
                            end
                            
                        end
                        
                        % book keeping & calculate focus index
                        stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
                        currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:),...
                            machine.meta.LUT_bxWidthminFWHM(2,:),...
                            pln.bixelWidth);
                        focusIx  =  ones(stf(i).numOfBixelsPerRay(j),1);
                        [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                            repmat(stf(i).ray(j).energy,length([machine.data]),1))));
                        
                        % get for each spot the focus index
                        for k = 1:stf(i).numOfBixelsPerRay(j)
                            focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                        end
                        
                        stf(i).ray(j).focusIx = focusIx';
                        
                    else % target not hit
                        stf(i).ray(j)               = [];
                        stf(i).numOfBixelsPerRay(j) = [];
                    end
                    
                elseif strcmp(stf(i).radiationMode,'photons')
                    
                    % book keeping for photons
                    stf(i).ray(j).energy = machine.data.energy;
                    
                else
                    error('Error generating stf struct: invalid radiation modality.');
                end
        end
        
        
        % store total number of rays for beam-i
        stf(i).numOfRays = size(stf(i).ray,2);
        
        % post processing for particle remove energy slices
        if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')
            
            % get minimum energy per field
            minEnergy = min([stf(i).ray.energy]);
            maxEnergy = max([stf(i).ray.energy]);
            
            % get corresponding peak position
            availableEnergies = [machine.data.energy];
            minPeakPos  = machine.data(minEnergy == availableEnergies).peakPos;
            maxPeakPos  = machine.data(maxEnergy == availableEnergies).peakPos;
            
            % find set of energyies with adequate spacing
            
            if strcmp(machine.meta.machine,'Generic')
                longitudinalSpotSpacing = 1.5; % enforce all entries to be used
            else
                longitudinalSpotSpacing = 3;   % default value for all other treatment machines
            end
            
            tolerance              = longitudinalSpotSpacing/10;
            availablePeakPos       = [machine.data.peakPos];
            
            useEnergyBool = availablePeakPos >= minPeakPos & availablePeakPos <= maxPeakPos;
            
            ixCurr = find(useEnergyBool,1,'first');
            ixRun  = ixCurr + 1;
            ixEnd  = find(useEnergyBool,1,'last');
            
            while ixRun <= ixEnd
                if abs(availablePeakPos(ixRun)-availablePeakPos(ixCurr)) < ...
                        longitudinalSpotSpacing - tolerance
                    useEnergyBool(ixRun) = 0;
                else
                    ixCurr = ixRun;
                end
  
                % book keeping & calculate focus index
                stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
                currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:)',...
                                             machine.meta.LUT_bxWidthminFWHM(2,:)',...
                                             pln.propStf.bixelWidth);
                focusIx  =  ones(stf(i).numOfBixelsPerRay(j),1);
                [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                                repmat(stf(i).ray(j).energy,length([machine.data]),1))));

                % get for each spot the focus index
                for k = 1:stf(i).numOfBixelsPerRay(j)                    
                    focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                end
                if isempty(stf(i).ray(j).energy)
                    stf(i).ray(j) = [];
                    stf(i).numOfBixelsPerRay(j) = [];
                    stf(i).numOfRays = stf(i).numOfRays - 1;
                end
            end
            
        end
        
        % save total number of bixels
        stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
        
    else
        %Determine which initialized beam the current beam belongs
        %to
        [~,stf(i).beamParentInitIndex] = min(abs(pln.propStf.initGantryAngles-pln.propStf.gantryAngles(i)));
        stf(i).beamParentGantryAngle = pln.propStf.initGantryAngles(stf(i).beamParentInitIndex);
        stf(i).beamParentIndex = find(pln.propStf.gantryAngles == stf(i).beamParentGantryAngle);
        
        %Indicate if this beam is to be included in optimization/initialization or not
        %All beams are still considered in dose calc for objective function
        stf(i).initializeBeam = any(pln.propStf.initGantryAngles==pln.propStf.gantryAngles(i));
        stf(i).optimizeBeam = any(pln.propStf.optGantryAngles==pln.propStf.gantryAngles(i));
        
        %Determine different angle borders
        %doseAngleBorders are the angular borders over which dose is
        %deposited
        
        if i == 1
            
            stf(i).doseAngleBorders = ([pln.propStf.gantryAngles(i) pln.propStf.gantryAngles(i+1)]+pln.propStf.gantryAngles(i))/2;
        elseif i == length(pln.propStf.gantryAngles)
            
            stf(i).doseAngleBorders = ([pln.propStf.gantryAngles(i-1) pln.propStf.gantryAngles(i)]+pln.propStf.gantryAngles(i))/2;
        else
            
            stf(i).doseAngleBorders = ([pln.propStf.gantryAngles(i-1) pln.propStf.gantryAngles(i+1)]+pln.propStf.gantryAngles(i))/2;
        end
        
        stf(i).doseAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).doseAngleBorders(1) stf(i).doseAngleBorders(2)-stf(i).gantryAngle];
        stf(i).doseAngleBordersDiff = sum(stf(i).doseAngleBorderCentreDiff);
        
        %Assign beam to its Parent, either as child (optimized) or subchild
        %(interpolated)
        if stf(i).optimizeBeam
            if ~isfield(stf(stf(i).beamParentIndex),'beamChildrenGantryAngles') || isempty(stf(stf(i).beamParentIndex).beamChildrenGantryAngles)
                stf(stf(i).beamParentIndex).numOfBeamChildren = 0;
                stf(stf(i).beamParentIndex).beamChildrenGantryAngles = nan(1000,1);
                stf(stf(i).beamParentIndex).beamChildrenIndex = nan(1000,1);
            end
            
            stf(stf(i).beamParentIndex).numOfBeamChildren = stf(stf(i).beamParentIndex).numOfBeamChildren+1;
            stf(stf(i).beamParentIndex).beamChildrenGantryAngles(stf(stf(i).beamParentIndex).numOfBeamChildren) = pln.propStf.gantryAngles(i);
            stf(stf(i).beamParentIndex).beamChildrenIndex(stf(stf(i).beamParentIndex).numOfBeamChildren) = i;
            
            %Determine different angle borders
            %optAngleBorders are the angular borders over which an optimized control point
            %has influence
            optIndex = find(pln.propStf.optGantryAngles == pln.propStf.gantryAngles(i));
            
            if optIndex == 1
                stf(i).optAngleBorders = ([pln.propStf.optGantryAngles(optIndex) pln.propStf.optGantryAngles(optIndex+1)]+pln.propStf.optGantryAngles(optIndex))/2;
                
                lastOptIndex = i;
                nextOptIndex = find(pln.propStf.gantryAngles == pln.propStf.optGantryAngles(optIndex+1));
                
                stf(i).lastOptIndex = i;
                stf(i).nextOptIndex = find(pln.propStf.gantryAngles == pln.propStf.optGantryAngles(optIndex+1));
            elseif optIndex == length(pln.propStf.optGantryAngles)
                stf(i).optAngleBorders = ([pln.propStf.optGantryAngles(optIndex-1) pln.propStf.optGantryAngles(optIndex)]+pln.propStf.optGantryAngles(optIndex))/2;
                
                stf(i).lastOptIndex = find(pln.propStf.gantryAngles == pln.propStf.optGantryAngles(optIndex-1));
                stf(i).nextOptIndex = i;
            else
                stf(i).optAngleBorders = ([pln.propStf.optGantryAngles(optIndex-1) pln.propStf.optGantryAngles(optIndex+1)]+pln.propStf.optGantryAngles(optIndex))/2;
                
                lastOptIndex = i;
                nextOptIndex = find(pln.propStf.gantryAngles == pln.propStf.optGantryAngles(optIndex+1));
                
                stf(i).lastOptIndex = find(pln.propStf.gantryAngles == pln.propStf.optGantryAngles(optIndex-1));
                stf(i).nextOptIndex = find(pln.propStf.gantryAngles == pln.propStf.optGantryAngles(optIndex+1));
            end
            stf(i).doseAngleOpt = ones(1,2);
            
            stf(i).optAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).optAngleBorders(1) stf(i).optAngleBorders(2)-stf(i).gantryAngle];
            stf(i).optAngleBordersDiff = sum(stf(i).optAngleBorderCentreDiff);
            
            %This is the factor that relates the total time in the
            %optimized arc sector to the total time in the current dose
            %sector
            stf(i).timeFacCurr =  stf(i).doseAngleBordersDiff./stf(i).optAngleBordersDiff;
            
            %These are the factors that relate the total time in the
            %optimized arc sector to the total time in the previous and
            %next dose sectors
            
        else
            if ~isfield(stf(stf(i).beamParentIndex),'beamSubChildrenGantryAngles') || isempty(stf(stf(i).beamParentIndex).beamSubChildrenGantryAngles)
                stf(stf(i).beamParentIndex).numOfBeamSubChildren = 0;
                stf(stf(i).beamParentIndex).beamSubChildrenGantryAngles = nan(1000,1);
                stf(stf(i).beamParentIndex).beamSubChildrenIndex = nan(1000,1);
            end
            
            stf(stf(i).beamParentIndex).numOfBeamSubChildren = stf(stf(i).beamParentIndex).numOfBeamSubChildren+1;
            stf(stf(i).beamParentIndex).beamSubChildrenGantryAngles(stf(stf(i).beamParentIndex).numOfBeamSubChildren) = pln.propStf.gantryAngles(i);
            stf(stf(i).beamParentIndex).beamSubChildrenIndex(stf(stf(i).beamParentIndex).numOfBeamSubChildren) = i;
            
            stf(i).fracFromLastOpt = (pln.propStf.gantryAngles(nextOptIndex)-pln.propStf.gantryAngles(i))./(pln.propStf.gantryAngles(nextOptIndex)-pln.propStf.gantryAngles(lastOptIndex));
            stf(i).lastOptIndex = lastOptIndex;
            stf(i).nextOptIndex = nextOptIndex;
        end
        
        
        if stf(i).initializeBeam
            %Determine different angle borders
            %initAngleBorders are the angular borders over which an optimized control point
            %has influence
            initIndex = find(pln.propStf.initGantryAngles == pln.propStf.gantryAngles(i));
            
            if initIndex == 1
                
                stf(i).initAngleBorders = [min(pln.propStf.initGantryAngles(initIndex),pln.propStf.gantryAngles(1)) (pln.propStf.initGantryAngles(initIndex+1)+pln.propStf.initGantryAngles(initIndex))/2];
            elseif initIndex == length(pln.propStf.initGantryAngles)
                
                stf(i).initAngleBorders = [(pln.propStf.initGantryAngles(initIndex-1)+pln.propStf.initGantryAngles(initIndex))/2 max(pln.propStf.initGantryAngles(initIndex),pln.propStf.gantryAngles(end))];
            else
                
                stf(i).initAngleBorders = ([pln.propStf.initGantryAngles(initIndex-1) pln.propStf.initGantryAngles(initIndex+1)]+pln.propStf.initGantryAngles(initIndex))/2;
            end
            stf(i).initAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).initAngleBorders(1) stf(i).initAngleBorders(2)-stf(i).gantryAngle];
            stf(i).initAngleBordersDiff = sum(stf(i).initAngleBorderCentreDiff);
        end
        
        %The following must be taken as the union of stf(:).FIELD and stf(:).FIELD:
        %ray.rayPos_bev
        %ray.targetPoint_bev
        %Then these are rotated to form the non-bev forms;
        %ray.rayCorners_SCD is also formed
        numOfRays = stf(i).numOfRays;
        rayPosBEV = reshape([stf(i).ray(:).rayPos_bev]',3,numOfRays)';
        targetPointBEV = reshape([stf(i).ray(:).targetPoint_bev]',3,numOfRays)';
        
        masterRayPosBEV = union(masterRayPosBEV,rayPosBEV,'rows');
        masterTargetPointBEV = union(masterTargetPointBEV,targetPointBEV,'rows');
    end
    
    
    
    % Show progress
    matRad_progress(i,length(pln.propStf.gantryAngles));
    
    %% visualization
    if visMode > 0
        
        clf;
        % first subplot: visualization in bev
        subplot(1,2,1)
        hold on
        
        % plot rotated target coordinates
        plot3(rot_coords(:,1),rot_coords(:,2),rot_coords(:,3),'r.')
        
        % surface rendering
        if visMode == 2
            
            % generate a 3D rectangular grid centered at isocenter in
            % voxel coordinates
            [X,Y,Z] = meshgrid((1:ct.cubeDim(2))-stf(i).isoCenter(1)/ct.resolution.x, ...
                (1:ct.cubeDim(1))-stf(i).isoCenter(2)/ct.resolution.y, ...
                (1:ct.cubeDim(3))-stf(i).isoCenter(3)/ct.resolution.z);
            
            % computes surface
            patSurfCube      = 0*ct.cube{1};
            idx              = [cst{:,4}];
            idx              = unique(vertcat(idx{:}));
            patSurfCube(idx) = 1;
            
            [f,v] = isosurface(X,Y,Z,patSurfCube,.5);
            
            % convert isosurface from voxel to [mm]
            v(:,1) = v(:,1)*ct.resolution.x;
            v(:,2) = v(:,2)*ct.resolution.y;
            v(:,3) = v(:,3)*ct.resolution.z;
            
            % rotate surface
            rotated_surface = v*rotMat_system_T;
            
            % surface rendering
            surface = patch('Faces',f,'Vertices',rotated_surface);
            set(surface,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',.4);
            lighting gouraud;
            
        end
        
        % plot projection matrix: coordinates at isocenter
        plot3(rayPos(:,1),rayPos(:,2),rayPos(:,3),'k.');
        
        % Plot matrix border of matrix at isocenter
        for j = 1:stf(i).numOfRays
            
            % Compute border for every bixels
            targetPoint_vox_X_1 = stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth;
            targetPoint_vox_Y_1 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_1 = stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth;
            
            targetPoint_vox_X_2 = stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth;
            targetPoint_vox_Y_2 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_2 = stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth;
            
            targetPoint_vox_X_3 = stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth;
            targetPoint_vox_Y_3 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_3 = stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth;
            
            targetPoint_vox_X_4 = stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth;
            targetPoint_vox_Y_4 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_4 = stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth;
            
            % plot
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_1],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_1],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_1],'g')
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_2],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_2],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_2],'g')
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_3],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_3],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_3],'g')
            plot3([stf(i).sourcePoint_bev(1) targetPoint_vox_X_4],[stf(i).sourcePoint_bev(2) targetPoint_vox_Y_4],[stf(i).sourcePoint_bev(3) targetPoint_vox_Z_4],'g')
            
        end
        
        % Plot properties
        daspect([1 1 1]);
        view(0,-90);
        xlabel 'X [mm]'
        ylabel 'Y [mm]'
        zlabel 'Z [mm]'
        title ('Beam''s eye view')
        axis([-300 300 -300 300 -300 300]);
        
        % second subplot: visualization in lps coordinate system
        subplot(1,2,2)
        
        % Plot target coordinates whitout any rotation
        plot3(coordsX,coordsY,coordsZ,'r.')
        hold on;
        
        % Rotated projection matrix at isocenter
        isocenter_plane_coor = rayPos*rotMat_vectors_T;
        
        % Plot isocenter plane
        plot3(isocenter_plane_coor(:,1),isocenter_plane_coor(:,2),isocenter_plane_coor(:,3),'y.');
        
        % Plot rotated bixels border.
        for j = 1:stf(i).numOfRays
            % Generate rotated projection target points.
            targetPoint_vox_1_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth]*rotMat_vectors_T;
            targetPoint_vox_2_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth]*rotMat_vectors_T;
            targetPoint_vox_3_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.propStf.bixelWidth]*rotMat_vectors_T;
            targetPoint_vox_4_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.propStf.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.propStf.bixelWidth]*rotMat_vectors_T;
            
            % Plot rotated target points.
            plot3([stf(i).sourcePoint(1) targetPoint_vox_1_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_1_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_1_rotated(:,3)],'g')
            plot3([stf(i).sourcePoint(1) targetPoint_vox_2_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_2_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_2_rotated(:,3)],'g')
            plot3([stf(i).sourcePoint(1) targetPoint_vox_3_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_3_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_3_rotated(:,3)],'g')
            plot3([stf(i).sourcePoint(1) targetPoint_vox_4_rotated(:,1)],[stf(i).sourcePoint(2) targetPoint_vox_4_rotated(:,2)],[stf(i).sourcePoint(3) targetPoint_vox_4_rotated(:,3)],'g')
        end
        
        % surface rendering
        if visMode == 2
            surface = patch('Faces',f,'Vertices',v);
            set(surface,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',.4);
            lighting gouraud;
        end
        
        % labels etc.
        daspect([1 1 1]);
        view(0,-90);
        xlabel 'X [mm]'
        ylabel 'Y [mm]'
        zlabel 'Z [mm]'
        title 'lps coordinate system'
        axis([-300 300 -300 300 -300 300]);
        %pause(1);
    end
    
    % include rangeshifter data if not yet available
    if strcmp(pln.radiationMode, 'protons') || strcmp(pln.radiationMode, 'carbon')
        for j = 1:stf(i).numOfRays
            for k = 1:numel(stf(i).ray(j).energy)
                stf(i).ray(j).rangeShifter(k).ID = 0;
                stf(i).ray(j).rangeShifter(k).eqThickness = 0;
                stf(i).ray(j).rangeShifter(k).sourceRashiDistance = 0;
            end
        end
    end
    
end

%% VMAT

if pln.propOpt.runVMAT
    numSSDErr = 0;
    
    %After all steering file information is completed, loop over
    %initialized gantry angles.  All children and subchildren of these angles should
    %have ray positions given by the union of their own ray positions and
    %the ray positions of the parent transformed to their gantry angle.
    %This is so that: (1) the target is still totally in the FOV of each
    %angle; and (2) the parent can give segments to the children during
    %initial segmentation and DAO.
    
    fprintf('matRad: Combining parent and child ray vectors in stf (VMAT)... ');
    
    
    for i = 1:length(pln.propStf.gantryAngles)
        if stf(i).optimizeBeam
            
            stf(i).timeFac = zeros(1,2);
            
            if i == 1
                stf(i).timeFac(1) = 0;
                stf(i).timeFac(2) = stf(i).optAngleBorderCentreDiff(2)/stf(i).optAngleBordersDiff;
            elseif i == length(pln.propStf.gantryAngles)
                stf(i).timeFac(1) = stf(i).optAngleBorderCentreDiff(1)/stf(i).optAngleBordersDiff;
                stf(i).timeFac(2) = 0;
            else
                stf(i).timeFac(1) = stf(i).optAngleBorderCentreDiff(1)/stf(i).optAngleBordersDiff;
                stf(i).timeFac(2) = stf(i).optAngleBorderCentreDiff(2)/stf(i).optAngleBordersDiff;
            end
            
            
            if stf(i).initializeBeam
                %remove NaNs from beamChildren and beamSubChildren
                if isfield(stf(i),'beamChildrenGantryAngles')
                    stf(i).beamChildrenGantryAngles(isnan(stf(i).beamChildrenGantryAngles)) = [];
                    stf(i).beamChildrenIndex(isnan(stf(i).beamChildrenIndex)) = [];
                else
                    stf(i).numOfBeamChildren = 0;
                end
                if isfield(stf(i),'beamSubChildrenGantryAngles')
                    stf(i).beamSubChildrenGantryAngles(isnan(stf(i).beamSubChildrenGantryAngles)) = [];
                    stf(i).beamSubChildrenIndex(isnan(stf(i).beamSubChildrenIndex)) = [];
                else
                    stf(i).numOfBeamSubChildren = 0;
                end
            end
        else
            
            % for time interpolation
            stf(i).timeFracFromLastOpt = (stf(stf(i).lastOptIndex).optAngleBorders(2)-stf(i).doseAngleBorders(1))./stf(i).doseAngleBordersDiff;
            stf(i).timeFracFromNextOpt = (stf(i).doseAngleBorders(2)-stf(stf(i).lastOptIndex).optAngleBorders(2))./stf(i).doseAngleBordersDiff;
            if stf(i).timeFracFromLastOpt > 1
                stf(i).timeFracFromLastOpt = 1;
            elseif stf(i).timeFracFromLastOpt < 0
                stf(i).timeFracFromLastOpt = 0;
            end
            if stf(i).timeFracFromNextOpt > 1
                stf(i).timeFracFromNextOpt = 1;
            elseif stf(i).timeFracFromNextOpt < 0
                stf(i).timeFracFromNextOpt = 0;
            end
        end
        
        currMasterRayPosBEV = masterRayPosBEV;
        currMasterTargetPointBEV = masterTargetPointBEV;
        
        currMasterRayPosBEV(isnan(currMasterRayPosBEV)) = [];
        currMasterTargetPointBEV(isnan(currMasterTargetPointBEV)) = [];
        currMasterRayPosBEV = reshape(currMasterRayPosBEV,[],3);
        currMasterTargetPointBEV = reshape(currMasterTargetPointBEV,[],3);
        
        stf(i).numOfRays = size(currMasterRayPosBEV,1);
        stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
        stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
        
        
        % source position in bev
        stf(i).sourcePoint_bev = [0 -SAD 0];
        
        % get (active) rotation matrix
        % transpose matrix because we are working with row vectors
        rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i)));
        
        
        stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;
        
        % Save ray and target position in lps system.
        for j = 1:stf(i).numOfRays
            stf(i).ray(j).rayPos_bev = currMasterRayPosBEV(j,:);
            stf(i).ray(j).targetPoint_bev = currMasterTargetPointBEV(j,:);
            
            stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMat_vectors_T;
            stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMat_vectors_T;
            if strcmp(pln.radiationMode,'photons')
                stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                    [currMasterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    currMasterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    currMasterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                    currMasterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMat_vectors_T;
            end
        end
        
        % loop over all rays to determine meta information for each ray
        stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
        
        for j = stf(i).numOfRays:-1:1
            
            % ray tracing necessary to determine depth of the target
            [~,l,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                ct.resolution, ...
                stf(i).sourcePoint, ...
                stf(i).ray(j).targetPoint, ...
                [ct.cube {voiTarget}]);
                
                % find appropriate energies for particles
                if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')
                    
                    % target hit
                    if sum(rho{2}) > 0
                        
                        % compute radiological depths
                        % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                        radDepths = cumsum(l .* rho{1});
                        
                        % find target entry & exit
                        diff_voi    = diff([rho{2}]);
                        targetEntry = radDepths(diff_voi == 1);
                        targetExit  = radDepths(diff_voi == -1);
                        
                        if numel(targetEntry) ~= numel(targetExit)
                            error('Inconsistency during ray tracing.');
                        end
                        
                        stf(i).ray(j).energy = [];
                        
                        % Save energies in stf struct
                        for k = 1:numel(targetEntry)
                            stf(i).ray(j).energy = [stf(i).ray(j).energy availableEnergies(availablePeakPos>=targetEntry(k)&availablePeakPos<=targetExit(k))];
                            % adjust spot spacing according to pln.bixelWidth when using HIT basedata
                            %DefaultLongitudialSpotSpacing = pln.bixelWidth;  % in [mm]
                            DefaultLongitudialSpotSpacing = 3;
                            if strcmp(pln.machine,'HIT') && length(stf(i).ray(j).energy)>2
                                Tolerance = 0.5;
                                hasNext = true;
                                CntEnergy =2;
                                while hasNext
                                    if abs(stf(i).ray(j).energy(CntEnergy)-stf(i).ray(j).energy(CntEnergy-1))<...
                                            DefaultLongitudialSpotSpacing-Tolerance
                                        stf(i).ray(j).energy(CntEnergy)=[];
                                    else
                                        CntEnergy = CntEnergy+1;
                                    end
                                    if CntEnergy == length(stf(i).ray(j).energy)
                                        hasNext = false;
                                    end
                                end
                            end
                            
                        end
                        
                        % book keeping & calculate focus index
                        stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
                        currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:),...
                            machine.meta.LUT_bxWidthminFWHM(2,:),...
                            pln.bixelWidth);
                        focusIx  =  ones(stf(i).numOfBixelsPerRay(j),1);
                        [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                            repmat(stf(i).ray(j).energy,length([machine.data]),1))));
                        
                        % get for each spot the focus index
                        for k = 1:stf(i).numOfBixelsPerRay(j)
                            focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                        end
                        
                        stf(i).ray(j).focusIx = focusIx';
                        
                    else % target not hit
                        stf(i).ray(j)               = [];
                        stf(i).numOfBixelsPerRay(j) = [];
                    end
                    
                elseif strcmp(stf(i).radiationMode,'photons')
                    
                    % book keeping for photons
                    stf(i).ray(j).energy = machine.data.energy;
                    
                else
                    error('Error generating stf struct: invalid radiation modality.');
                end
        end
        
        matRad_progress(i,length(pln.propStf.gantryAngles));
    end
    
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
    
    fprintf('\n\nNumber of SSD errors: %d out of %d.\n\n',numSSDErr,sum([stf(:).numOfRays]));
end

end

