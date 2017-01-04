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

if numel(pln.gantryAngles) ~= numel(pln.couchAngles)
    error('Inconsistent number of gantry and couch angles.');
end

if pln.bixelWidth < 0 || ~isfinite(pln.bixelWidth)
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

% density threshold used for SSD calculation
DensityThresholdSSD = 0.05;

% generate voi cube for targets
voiTarget    = zeros(ct.cubeDim);
voiTarget(V) = 1;
    
% add margin
addmarginBool = 1;
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

% Correct for iso center position. Whit this correction Isocenter is
% (0,0,0) [mm]
coordsX = coordsX_vox*ct.resolution.x - pln.isoCenter(1);
coordsY = coordsY_vox*ct.resolution.y - pln.isoCenter(2);
coordsZ = coordsZ_vox*ct.resolution.z - pln.isoCenter(3);

% Define steering file like struct. Prellocating for speed.
stf = struct;

if pln.VMAT
    %Initialize master ray positions and target points with NaNs, to be
    %deleted later.  These arrays are the unions of the corresponding
    %arrays per gantry angle.  In order to do VMAT, it is easier to have
    %the same MLC range and dij calculation for every possible beam/gantry
    %angle.
    masterRayPosBEV = nan(1,3);
    masterTargetPointBEV = nan(1,3);
    
    stf(1).defaultGantryRot = pln.defaultGantryRot;
end
    

% loop over all angles
for i = 1:length(pln.gantryAngles)
    
    % Save meta information for treatment plan
    stf(i).gantryAngle   = pln.gantryAngles(i);
    stf(i).couchAngle    = pln.couchAngles(i);
    stf(i).bixelWidth    = pln.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    stf(i).SAD           = SAD;
    stf(i).isoCenter     = pln.isoCenter;
    
    % gantry and couch roation matrices according to IEC 61217 standard
    % instead of moving the beam around the patient, we perform an inverse
    % rotation of the patient, i.e. we consider a beam's eye view
    % coordinate system; use transpose matrices because we are working with
    % row vectors
    
    % Rotation around Z axis (gantry)
    inv_rotMx_XY_T = [ cosd(-pln.gantryAngles(i)) sind(-pln.gantryAngles(i)) 0;
                      -sind(-pln.gantryAngles(i)) cosd(-pln.gantryAngles(i)) 0;
                                                0                          0 1];
    
    % Rotation around Y axis (Couch movement)
    inv_rotMx_XZ_T = [cosd(-pln.couchAngles(i)) 0 -sind(-pln.couchAngles(i));
                                              0 1                          0;
                      sind(-pln.couchAngles(i)) 0  cosd(-pln.couchAngles(i))];
    
    % rotate target coordinates (1st couch around Y axis, 2nd gantry around
    % z axis); matrix multiplication not cummutative
    rot_coords = [coordsX coordsY coordsZ]*inv_rotMx_XZ_T*inv_rotMx_XY_T;
    
    % project x and z coordinates to isocenter
    coordsAtIsoCenterPlane(:,1) = (rot_coords(:,1)*SAD)./(SAD + rot_coords(:,2));
    coordsAtIsoCenterPlane(:,2) = (rot_coords(:,3)*SAD)./(SAD + rot_coords(:,2));
    
    % Take unique rows values for beamlets positions. Calculate position of
    % central ray for every bixel    
    rayPos = unique(pln.bixelWidth*round([            coordsAtIsoCenterPlane(:,1) ... 
                                          zeros(size(coordsAtIsoCenterPlane,1),1) ...
                                                      coordsAtIsoCenterPlane(:,2)]/pln.bixelWidth),'rows');
                                                  
    % pad ray position array if resolution of target voxel grid not sufficient
    maxCtResolution = max([ct.resolution.x ct.resolution.y ct.resolution.z]);
    if pln.bixelWidth < maxCtResolution
        origRayPos = rayPos;
        for j = -floor(maxCtResolution/pln.bixelWidth):floor(maxCtResolution/pln.bixelWidth)
            for k = -floor(maxCtResolution/pln.bixelWidth):floor(maxCtResolution/pln.bixelWidth)
                if abs(j)+abs(k)==0
                    continue;
                end                
                rayPos = [rayPos; origRayPos(:,1)+j*pln.bixelWidth origRayPos(:,2) origRayPos(:,3)+k*pln.bixelWidth];
            end
        end
     end

     % remove spaces within rows of bixels for DAO
     if pln.runDAO
         % create single x,y,z vectors
         x = rayPos(:,1);
         y = rayPos(:,2);
         z = rayPos(:,3);
         uniZ = unique(z);
         for j = 1:numel(uniZ)
             x_loc = x(z == uniZ(j));
             x_min = min(x_loc);
             x_max = max(x_loc);
             x = [x; [x_min:pln.bixelWidth:x_max]'];
             y = [y; zeros((x_max-x_min)/pln.bixelWidth+1,1)];
             z = [z; uniZ(j)*ones((x_max-x_min)/pln.bixelWidth+1,1)];             
         end
         
         rayPos = [x,y,z];
     end
    
    % remove double rays
    rayPos = unique(rayPos,'rows');
    
    % Save the number of rays
    stf(i).numOfRays = size(rayPos,1);
    
    % Save ray and target position in beam eye´s view (bev)
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos_bev = rayPos(j,:);
        stf(i).ray(j).targetPoint_bev = [2*stf(i).ray(j).rayPos_bev(1) ...
                                                               SAD ...
                                         2*stf(i).ray(j).rayPos_bev(3)];
    end
    
    % source position in bev
    stf(i).sourcePoint_bev = [0 -SAD 0];
    
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient; use transpose matrices because we are
    % working with row vectors
    
    % Rotation around Z axis (gantry)
    rotMx_XY_T = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
                  -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
                                           0                         0 1];
    
    % Rotation around Y axis (couch)
    rotMx_XZ_T = [cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
                                         0 1                        0;
                  sind(pln.couchAngles(i)) 0  cosd(pln.couchAngles(i))];
    
    % Rotated Source point (1st gantry, 2nd couch)
    stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMx_XY_T*rotMx_XZ_T;
    
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMx_XY_T*rotMx_XZ_T;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMx_XY_T*rotMx_XZ_T;
        stf(i).ray(j).SSD         = NaN;
        if strcmp(pln.radiationMode,'photons') 
            stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                                                             [rayPos(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                                                              rayPos(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                                                              rayPos(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                                                              rayPos(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMx_XY_T*rotMx_XZ_T;
        end
    end
    
    % loop over all rays to determine meta information for each ray    
    stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
    
    if ~pln.VMAT
        %If it is VMAT, we will do this later
        for j = stf(i).numOfRays:-1:1
            
            % ray tracing necessary to determine depth of the target
            [alpha,l,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                ct.resolution, ...
                stf(i).sourcePoint, ...
                stf(i).ray(j).targetPoint, ...
                [ct.cube {voiTarget}]);
            
            ixSSD = find(rho{1} > DensityThresholdSSD,1,'first');
            
            if isempty(ixSSD)== 1
                warning('Surface for SSD calculation starts directly in first voxel of CT\n');
            end
            
            % calculate SSD
            stf(i).ray(j).SSD = 2 * stf(i).SAD * alpha(ixSSD);
            
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
    end
    
    % store total number of rays for beam-i
    stf(i).numOfRays = size(stf(i).ray,2);
     
    % save total number of bixels
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    %     figure,
    %     for jj = 1:length(stf.ray)
    %        plot(stf.ray(jj).rayPos_bev(1),stf.ray(jj).rayPos_bev(3),'rx'),hold on
    %     end
    
    
    if pln.VMAT
        %include min/max gantry rotation speed in first steering file
        stf(1).gantryRotCst = pln.gantryRotCst;
        %Indicate if this beam is to be included in optimization/initialization or not
        %All beams are still considered in dose calc for objective function
        numInitGantryAngles = length(pln.initGantryAngles);
        
        if any(pln.initGantryAngles==pln.gantryAngles(i)) %Assume all initialized gantry angles are included in DAO
            stf(i).initializeBeam = true;
            stf(i).optimizeBeam = true;
            stf(i).beamChildrenIndex = zeros(10000,1);
            stf(i).beamChildrenGantryAngles = zeros(10000,1);
            stf(i).numOfBeamChildren = 0;
            stf(i).beamSubChildrenIndex = zeros(10000,1);
            stf(i).beamSubChildrenGantryAngles = zeros(10000,1);
            stf(i).numOfBeamSubChildren = 0;
            
            initIndex = find(pln.initGantryAngles == pln.gantryAngles(i));
            
            currInitGantryAngle = pln.initGantryAngles(initIndex);
            if currInitGantryAngle == 0
                nextInitGantryAngle = pln.initGantryAngles(mod0(initIndex+1,numInitGantryAngles));
                prevInitGantryAngle = -0.01;
            elseif currInitGantryAngle == 360
                nextInitGantryAngle = 360.01;
                prevInitGantryAngle = pln.initGantryAngles(mod0(initIndex-1,numInitGantryAngles));
            else
                nextInitGantryAngle = pln.initGantryAngles(mod0(initIndex+1,numInitGantryAngles));
                prevInitGantryAngle = pln.initGantryAngles(mod0(initIndex-1,numInitGantryAngles));
            end

            closerToCurrGantryAngle = 1;
            itera = 0;
            iterb = 0;
            while closerToCurrGantryAngle
                %Perform search for optimized (and interpolated) gantry angles close to
                %initialized angles (in the sense of sequencing).
                j = mod0(i+ceil(itera/2)*(-1)^iterb,length(pln.gantryAngles));
                
                testGantryAngle = pln.gantryAngles(j);
                
                currDiff = abs(testGantryAngle-currInitGantryAngle);
                nextDiff = abs(testGantryAngle-nextInitGantryAngle);
                prevDiff = abs(testGantryAngle-prevInitGantryAngle);
                
                if itera == iterb && (currDiff < nextDiff && currDiff < prevDiff)
                    %Alternate between angles smaller than and greater than
                    %current initialized gantry angle to search for any that are closer
                    %to this one than the next
                    if any(pln.optGantryAngles == pln.gantryAngles(j))
                        stf(i).numOfBeamChildren = stf(i).numOfBeamChildren+1;
                        stf(i).beamChildrenIndex(stf(i).numOfBeamChildren) = j;
                        stf(i).beamChildrenGantryAngles(stf(i).numOfBeamChildren) = testGantryAngle;
                    else
                        stf(i).numOfBeamSubChildren = stf(i).numOfBeamSubChildren+1;
                        stf(i).beamSubChildrenIndex(stf(i).numOfBeamSubChildren) = j;
                        stf(i).beamSubChildrenGantryAngles(stf(i).numOfBeamSubChildren) = testGantryAngle;
                    end
                    
                    itera = itera+1;
                    iterb = iterb+1;
                elseif  itera == iterb && (currDiff == nextDiff || currDiff == prevDiff)
                    %Include the last one as a subchild
                    if any(pln.optGantryAngles == pln.gantryAngles(j))
                        stf(i).numOfBeamChildren = stf(i).numOfBeamChildren+1;
                        stf(i).beamChildrenIndex(stf(i).numOfBeamChildren) = j;
                        stf(i).beamChildrenGantryAngles(stf(i).numOfBeamChildren) = testGantryAngle;
                    else
                        stf(i).numOfBeamSubChildren = stf(i).numOfBeamSubChildren+1;
                        stf(i).beamSubChildrenIndex(stf(i).numOfBeamSubChildren) = j;
                        stf(i).beamSubChildrenGantryAngles(stf(i).numOfBeamSubChildren) = testGantryAngle;
                    end
                   
                    iterb = iterb+1;
                    border1 = testGantryAngle;
                    borderInd1 = j;
                elseif  itera == iterb && (currDiff > nextDiff || currDiff > prevDiff)
                    %Now only go in one direction by locking iterb (to the next value).  Do not assume angles
                    %are equally spaced out on both sides.
                    itera = itera+2*(1-mod(iterb,2));
                    % = itera+2 if iterb is even, itera if iterb is odd
                    % prevents double-counting certain j's if iterb is even
                    iterb = iterb+1;
                    border1 = testGantryAngle;
                    borderInd1 = j;
                elseif itera ~= iterb && (currDiff < nextDiff && currDiff < prevDiff)
                    if any(pln.optGantryAngles == pln.gantryAngles(j))
                        stf(i).numOfBeamChildren = stf(i).numOfBeamChildren+1;
                        stf(i).beamChildrenIndex(stf(i).numOfBeamChildren) = j;
                        stf(i).beamChildrenGantryAngles(stf(i).numOfBeamChildren) = testGantryAngle;
                    else
                        stf(i).numOfBeamSubChildren = stf(i).numOfBeamSubChildren+1;
                        stf(i).beamSubChildrenIndex(stf(i).numOfBeamSubChildren) = j;
                        stf(i).beamSubChildrenGantryAngles(stf(i).numOfBeamSubChildren) = testGantryAngle;
                    end
                    itera = itera+2; %skip forward 2 because of the ceil(itera/2)
                elseif itera ~= iterb && (currDiff >= nextDiff || currDiff >= prevDiff)
                    %Once we have exhausted close beam angles in both
                    %directions, finish search.
                    %{
                    if any(pln.optGantryAngles == pln.gantryAngles(j))
                        stf(i).numOfBeamChildren = stf(i).numOfBeamChildren+1;
                        stf(i).beamChildrenIndex(stf(i).numOfBeamChildren) = j;
                        stf(i).beamChildrenGantryAngles(stf(i).numOfBeamChildren) = testGantryAngle;
                    else
                        stf(i).numOfBeamSubChildren = stf(i).numOfBeamSubChildren+1;
                        stf(i).beamSubChildrenIndex(stf(i).numOfBeamSubChildren) = j;
                        stf(i).beamSubChildrenGantryAngles(stf(i).numOfBeamSubChildren) = testGantryAngle;
                    end
                    %}
                    closerToCurrGantryAngle = 0;
                    border2 = testGantryAngle;
                    borderInd2 = j;
                end
                
            end
            stf(i).beamChildrenGantryAngles(stf(i).beamChildrenIndex==0) = [];
            stf(i).beamChildrenIndex(stf(i).beamChildrenIndex==0) = [];
            [stf(i).beamChildrenGantryAngles,ind] = sortangles(stf(i).beamChildrenGantryAngles);
            stf(i).beamChildrenGantryAngles = unique(stf(i).beamChildrenGantryAngles,'stable');
            stf(i).beamChildrenIndex = stf(i).beamChildrenIndex(ind);
            stf(i).beamChildrenIndex = unique(stf(i).beamChildrenIndex,'stable');
            
            stf(i).beamSubChildrenGantryAngles(stf(i).beamSubChildrenIndex==0) = [];
            stf(i).beamSubChildrenIndex(stf(i).beamSubChildrenIndex==0) = [];
            [stf(i).beamSubChildrenGantryAngles,ind] = sortangles(stf(i).beamSubChildrenGantryAngles);
            stf(i).beamSubChildrenGantryAngles = unique(stf(i).beamSubChildrenGantryAngles,'stable');
            stf(i).beamSubChildrenIndex = stf(i).beamSubChildrenIndex(ind);
            stf(i).beamSubChildrenIndex = unique(stf(i).beamSubChildrenIndex,'stable');
            
            %borderAngles specifies the first and last angles in the arc
            %sector (usually subChildren, except for first and last)
            border1 = min([min(stf(i).beamSubChildrenGantryAngles) min(stf(i).beamChildrenGantryAngles)]);
            border2 = max([max(stf(i).beamSubChildrenGantryAngles) max(stf(i).beamChildrenGantryAngles)]);
            borderInd1 = min([min(stf(i).beamSubChildrenIndex) min(stf(i).beamChildrenIndex)]);
            borderInd2 = max([max(stf(i).beamSubChildrenIndex) max(stf(i).beamChildrenIndex)]);
            %stf(i).borderAngles = [border1 border2];
            stf(i).borderAnglesIndex = [borderInd1 borderInd2];
            
            %[stf(i).borderAngles,ind] = sortangles([border1 border2]);
            %borderInds = [borderInd1 borderInd2];
            %stf(i).borderAnglesIndex = borderInds(ind);
            
            optGantryAngles0to360 = pln.optGantryAngles;
            %optGantryAngles0to360(optGantryAngles0to360 == 0) = 360;
            for j = 1:numel(stf(i).beamSubChildrenIndex)
                stf(stf(i).beamSubChildrenIndex(j)).lastOptAngle = max(pln.optGantryAngles(0 <= stf(i).beamSubChildrenGantryAngles(j)-pln.optGantryAngles));
                stf(stf(i).beamSubChildrenIndex(j)).lastOptInd = find(pln.gantryAngles == stf(stf(i).beamSubChildrenIndex(j)).lastOptAngle);
                
                stf(stf(i).beamSubChildrenIndex(j)).nextOptAngle = min(pln.optGantryAngles(0 <= pln.optGantryAngles-stf(i).beamSubChildrenGantryAngles(j)));
                if isempty(stf(stf(i).beamSubChildrenIndex(j)).nextOptAngle)
                    stf(stf(i).beamSubChildrenIndex(j)).nextOptAngle = min(pln.optGantryAngles(0 <= optGantryAngles0to360-stf(i).beamSubChildrenGantryAngles(j)));
                end
                stf(stf(i).beamSubChildrenIndex(j)).nextOptInd = find(pln.gantryAngles == stf(stf(i).beamSubChildrenIndex(j)).nextOptAngle);
            end
           
            %Don't think commented part is necessary anymore since I
            %changed the first angle from 0->nonzero
            
            if stf(i).gantryAngle == pln.initGantryAngles(1)
                %                stf(i).borderAngles(1) = 2*stf(i).gantryAngle-stf(i).borderAngles(2);
                %                stf(i).borderAnglesIndex(1) = find(pln.gantryAngles == mod(stf(i).borderAngles(1),360));
                stf(i).borderAngles = [0 mean([currInitGantryAngle nextInitGantryAngle])];
                stf(i).lastBorderAngle = 0;
            elseif stf(i).gantryAngle == pln.initGantryAngles(end)
                %                stf(i).borderAngles(2) = 2*stf(i).gantryAngle-stf(i).borderAngles(1);
                %                stf(i).borderAnglesIndex(2) = find(pln.gantryAngles == mod0(stf(i).borderAngles(2),360));
                %stf(lastInitInd).nextBorderAngle = stf(i).borderAngles(1);
                %stf(i).lastBorderAngle = lastBorderAngle;
                %stf(i).nextBorderAngle = 360+stf(i).borderAngles(2);
                stf(i).borderAngles = [mean([prevInitGantryAngle currInitGantryAngle]) 360];
                stf(i).lastBorderAngle = lastBorderAngle;
                stf(i).nextBorderAngle = 360;
            else
                if exist('lastInitInd','var')
                    stf(i).borderAngles = [mean([prevInitGantryAngle currInitGantryAngle]) mean([currInitGantryAngle nextInitGantryAngle])];
                    stf(lastInitInd).nextBorderAngle = stf(i).borderAngles(1);
                    stf(i).lastBorderAngle = lastBorderAngle;
                end
            end
            lastInitInd = i;
            lastBorderAngle = stf(i).borderAngles(2);
            
            optIndex = find(pln.optGantryAngles == pln.gantryAngles(i));
            if optIndex == numel(pln.optGantryAngles)
                stf(i).nextOptAngleDiff = pln.optGantryAngles(optIndex)-pln.optGantryAngles(optIndex-1); %this is prev not next
            else
                stf(i).nextOptAngleDiff = pln.optGantryAngles(optIndex+1)-pln.optGantryAngles(optIndex);
            end
            if i == numel(pln.gantryAngles)
                stf(i).nextAngleDiff = pln.gantryAngles(i)-pln.gantryAngles(i-1); %this is prev not next
            else
                stf(i).nextAngleDiff = pln.gantryAngles(i+1)-pln.gantryAngles(i);
            end
            
        elseif any(pln.optGantryAngles==pln.gantryAngles(i)) && ~any(pln.initGantryAngles==pln.gantryAngles(i))
            %Some are not initialized, they are given apertures from the initialized angles, yet they are still independently optimized in DAO
            stf(i).initializeBeam = false;
            stf(i).optimizeBeam = true;
            optIndex = find(pln.optGantryAngles == pln.gantryAngles(i));
            if optIndex == numel(pln.optGantryAngles)
                stf(i).nextOptAngleDiff = pln.optGantryAngles(optIndex)-pln.optGantryAngles(optIndex-1); %this is prev not next
            else
                stf(i).nextOptAngleDiff = pln.optGantryAngles(optIndex+1)-pln.optGantryAngles(optIndex);
            end
            if i == numel(pln.gantryAngles)
                stf(i).nextAngleDiff = pln.gantryAngles(i)-pln.gantryAngles(i-1); %this is prev not next
            else
                stf(i).nextAngleDiff = pln.gantryAngles(i+1)-pln.gantryAngles(i);
            end
        else
            %These are purely interpolated in the DAO step, but dose is
            %calculated at these angles to improve dose calc accuracy.
            stf(i).initializeBeam = false;
            stf(i).optimizeBeam = false;
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
    matRad_progress(i,length(pln.gantryAngles));

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
            rotated_surface = v*inv_rotMx_XZ_T*inv_rotMx_XY_T;
            
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
            targetPoint_vox_X_1 = stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth;
            targetPoint_vox_Y_1 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_1 = stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth;
            
            targetPoint_vox_X_2 = stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth;
            targetPoint_vox_Y_2 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_2 = stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth;
            
            targetPoint_vox_X_3 = stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth;
            targetPoint_vox_Y_3 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_3 = stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth;
            
            targetPoint_vox_X_4 = stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth;
            targetPoint_vox_Y_4 = stf(i).ray(j).targetPoint_bev(:,2);
            targetPoint_vox_Z_4 = stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth;
            
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
        isocenter_plane_coor = rayPos*rotMx_XY_T*rotMx_XZ_T;
        
        % Plot isocenter plane
        plot3(isocenter_plane_coor(:,1),isocenter_plane_coor(:,2),isocenter_plane_coor(:,3),'y.');
        
        % Plot rotated bixels border.
        for j = 1:stf(i).numOfRays
            % Generate rotated projection target points.
            targetPoint_vox_1_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth]*rotMx_XY_T*rotMx_XZ_T;
            targetPoint_vox_2_rotated = [stf(i).ray(j).targetPoint_bev(:,1) + pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth]*rotMx_XY_T*rotMx_XZ_T;
            targetPoint_vox_3_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) - pln.bixelWidth]*rotMx_XY_T*rotMx_XZ_T;
            targetPoint_vox_4_rotated = [stf(i).ray(j).targetPoint_bev(:,1) - pln.bixelWidth,stf(i).ray(j).targetPoint_bev(:,2),stf(i).ray(j).targetPoint_bev(:,3) + pln.bixelWidth]*rotMx_XY_T*rotMx_XZ_T;
            
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
        
end    

%% VMAT
if pln.VMAT
    %After all steering file information is completed, loop over
    %initialized gantry angles.  All children and subchildren of these angles should
    %have ray positions given by the union of their own ray positions and
    %the ray positions of the parent transformed to their gantry angle.
    %This is so that: (1) the target is still totally in the FOV of each
    %angle; and (2) the parent can give segments to the children during
    %initial segmentation and DAO.
    
    fprintf('matRad: Combining parent and child ray vectors in stf (VMAT)... ');
    
    %Remove NaNs
    masterRayPosBEV(isnan(masterRayPosBEV)) = [];
    masterTargetPointBEV(isnan(masterTargetPointBEV)) = [];
    masterRayPosBEV = reshape(masterRayPosBEV,[],3);
    masterTargetPointBEV = reshape(masterTargetPointBEV,[],3);
    
    for i = 1:length(pln.gantryAngles)
        %This is a new rotation matrix, which is meant to
        %carry over the rayPos and
        %targetPoint from the optimized gantry angle to the
        %children gantry angles.
        stf(i).numOfRays = size(masterRayPosBEV,1);
        stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
        stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
        
        
        for j = 1:stf(i).numOfRays
            %rotate rayPos and targetPoint from BEV to patient coordinate
            %system
            rotMx_XY_T = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
                -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
                0                         0                         1];
            rotMx_XZ_T = [cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
                0                        1  0;
                sind(pln.couchAngles(i)) 0  cosd(pln.couchAngles(i))];
            
            stf(i).ray(j).rayPos_bev = masterRayPosBEV(j,:);
            stf(i).ray(j).targetPoint_bev = masterTargetPointBEV(j,:);
            
            
            stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMx_XY_T*rotMx_XZ_T;
            stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMx_XY_T*rotMx_XZ_T;
            if strcmp(pln.radiationMode,'photons')
                stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                    [masterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    masterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    masterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                    masterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMx_XY_T*rotMx_XZ_T;
            end
            
            % ray tracing necessary to determine depth of the target
            [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                ct.resolution, ...
                stf(i).sourcePoint, ...
                stf(i).ray(j).targetPoint, ...
                [ct.cube {voiTarget}]);
            
            ixSSD = find(rho{1} > DensityThresholdSSD,1,'first');
            
            if isempty(ixSSD)== 1
                warning('Surface for SSD calculation starts directly in first voxel of CT\n');
            end
            
            % calculate SSD
            stf(i).ray(j).SSD = 2 * stf(i).SAD * alpha(ixSSD);
            
            % book keeping for photons
            stf(i).ray(j).energy = machine.data.energy;
            
        end
        matRad_progress(i,length(pln.gantryAngles));
    end
    
    
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
    %{
    for parent = 1:length(pln.gantryAngles)
        if stf(parent).initializeBeam
            parentNumOfRays = stf(parent).numOfRays;
            parentRayPosBEV = reshape([stf(parent).ray(:).rayPos_bev]',3,parentNumOfRays)';
            parentTargetPointBEV = reshape([stf(parent).ray(:).targetPoint_bev]',3,parentNumOfRays)';
            
            for child = 1:stf(parent).numOfBeamChildren
                %The following must be taken as the union of stf(child).FIELD and stf(parent).FIELD:
                %ray.rayPos_bev
                %ray.targetPoint_bev
                %Then these are rotated to form the non-bev forms;
                %ray.rayCorners_SCD is also formed
                childIndex = stf(parent).beamChildrenIndex(child);
                childNumOfRays = stf(childIndex).numOfRays;
                
                childRayPosBEV = reshape([stf(childIndex).ray(:).rayPos_bev]',3,childNumOfRays)';
                childTargetPointBEV = reshape([stf(childIndex).ray(:).targetPoint_bev]',3,childNumOfRays)';
                
                childRayPosBEV_new = union(parentRayPosBEV,childRayPosBEV,'rows');
                childRayTargetPointBEV_new = union(parentTargetPointBEV,childTargetPointBEV,'rows');
                stf(childIndex).numOfRays = size(childRayPosBEV_new,1);
                
                %This is a new rotation matrix, which is meant to
                %carry over the rayPos and
                %targetPoint from the optimized gantry angle to the
                %children gantry angles.
                rotMx_XY_T = [ cosd(pln.gantryAngles(childIndex)) sind(pln.gantryAngles(childIndex)) 0;
                              -sind(pln.gantryAngles(childIndex)) cosd(pln.gantryAngles(childIndex)) 0;
                               0                                  0                                  1];
                rotMx_XZ_T = [cosd(pln.couchAngles(childIndex)) 0 -sind(pln.couchAngles(childIndex));
                              0                                 1  0;
                              sind(pln.couchAngles(childIndex)) 0  cosd(pln.couchAngles(childIndex))];
                
                for j = 1:stf(childIndex).numOfRays
                    stf(childIndex).ray(j).rayPos_bev = childRayPosBEV_new(j,:);
                    stf(childIndex).ray(j).targetPoint_bev = childRayTargetPointBEV_new(j,:);
                    
                    stf(childIndex).ray(j).rayPos      = stf(childIndex).ray(j).rayPos_bev*rotMx_XY_T*rotMx_XZ_T;
                    stf(childIndex).ray(j).targetPoint = stf(childIndex).ray(j).targetPoint_bev*rotMx_XY_T*rotMx_XZ_T;
                    if strcmp(pln.radiationMode,'photons')
                        stf(childIndex).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                            [childRayPosBEV_new(j,:) + [+stf(childIndex).bixelWidth/2,0,+stf(childIndex).bixelWidth/2];...
                            childRayPosBEV_new(j,:) + [-stf(childIndex).bixelWidth/2,0,+stf(childIndex).bixelWidth/2];...
                            childRayPosBEV_new(j,:) + [-stf(childIndex).bixelWidth/2,0,-stf(childIndex).bixelWidth/2];...
                            childRayPosBEV_new(j,:) + [+stf(childIndex).bixelWidth/2,0,-stf(childIndex).bixelWidth/2]])*rotMx_XY_T*rotMx_XZ_T;
                    end
                end
            end
            
            for subchild = 1:stf(parent).numOfBeamSubChildren
                %Do the exact same thing for the subchildren (interpolated
                %gantry angles).
                
                %The following must be taken as the union of stf(subchild).FIELD and stf(parent).FIELD:
                %ray.rayPos_bev
                %ray.targetPoint_bev
                %Then these are rotated to form the non-bev forms;
                %ray.rayCorners_SCD is also formed
                subChildIndex = stf(parent).beamSubChildrenIndex(subchild);
                subchildNumOfRays = stf(subChildIndex).numOfRays;
                
                subchildRayPosBEV = reshape([stf(subChildIndex).ray(:).rayPos_bev]',3,subchildNumOfRays)';
                subchildTargetPointBEV = reshape([stf(subChildIndex).ray(:).targetPoint_bev]',3,subchildNumOfRays)';
                
                subchildRayPosBEV_new = union(parentRayPosBEV,subchildRayPosBEV,'rows');
                subchildRayTargetPointBEV_new = union(parentTargetPointBEV,subchildTargetPointBEV,'rows');
                stf(subChildIndex).numOfRays = size(subchildRayPosBEV_new,1);
                
                %This is a new rotation matrix, which is meant to
                %carry over the rayPos and
                %targetPoint from the optimized gantry angle to the
                %subchildren gantry angles.
                rotMx_XY_T = [ cosd(pln.gantryAngles(subChildIndex)) sind(pln.gantryAngles(subChildIndex)) 0;
                              -sind(pln.gantryAngles(subChildIndex)) cosd(pln.gantryAngles(subChildIndex)) 0;
                               0                                     0                                     1];
                rotMx_XZ_T = [cosd(pln.couchAngles(subChildIndex)) 0 -sind(pln.couchAngles(subChildIndex));
                              0                                    1  0;
                              sind(pln.couchAngles(subChildIndex)) 0  cosd(pln.couchAngles(subChildIndex))];
                
                for j = 1:stf(subChildIndex).numOfRays
                    stf(subChildIndex).ray(j).rayPos_bev = subchildRayPosBEV_new(j,:);
                    stf(subChildIndex).ray(j).targetPoint_bev = subchildRayTargetPointBEV_new(j,:);
                    
                    stf(subChildIndex).ray(j).rayPos      = stf(subChildIndex).ray(j).rayPos_bev*rotMx_XY_T*rotMx_XZ_T;
                    stf(subChildIndex).ray(j).targetPoint = stf(subChildIndex).ray(j).targetPoint_bev*rotMx_XY_T*rotMx_XZ_T;
                    if strcmp(pln.radiationMode,'photons')
                        stf(subChildIndex).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                            [subchildRayPosBEV_new(j,:) + [+stf(subChildIndex).bixelWidth/2,0,+stf(subChildIndex).bixelWidth/2];...
                            subchildRayPosBEV_new(j,:) + [-stf(subChildIndex).bixelWidth/2,0,+stf(subChildIndex).bixelWidth/2];...
                            subchildRayPosBEV_new(j,:) + [-stf(subChildIndex).bixelWidth/2,0,-stf(subChildIndex).bixelWidth/2];...
                            subchildRayPosBEV_new(j,:) + [+stf(subChildIndex).bixelWidth/2,0,-stf(subChildIndex).bixelWidth/2]])*rotMx_XY_T*rotMx_XZ_T;
                    end
                end
            end
        end
        
        %Show progress
        matRad_progress(parent,length(pln.gantryAngles));
    end
    %}
end


end



function b = mod0(a,m)
%Modified version of MATLAB mod function, returns b=m where regular version
%would return 0 (i.e., a is 0 a multiple of m)

b = mod(a,m);
if b == 0
    b = m;
end

end

function [y,i] = sortangles(x)
%Modified version of MATLAB sort function, returns sorted angles,
%considering wraparound at \theta = 0, 360.

[y,i] = sort(x);
%if (max(y(:))-min(y(:))) >= 180
%    oldInd180 = find(y<=180,1,'last');
%    ln = length(y);
 %   newInd0 = ln-oldInd180+1;
 %   
 %   temp = zeros(size(y));
  %%  temp(1:(newInd0-1)) = y((oldInd180+1):ln);
  %  temp(newInd0:ln) = y(1:oldInd180);
  %  y = temp;
    
  %  tempi = zeros(size(y));
  %  tempi(1:(newInd0-1)) = i((oldInd180+1):ln);
  %%  tempi(newInd0:ln) = i(1:oldInd180);
   % i = tempi;
    
%end

end
