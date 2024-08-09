function stf = matRad_generateStf(ct,cst,pln,visMode)
% matRad steering information generation
%
% call
%   stf = matRad_generateStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting
%               this value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispInfo('matRad: Generating stf struct... ');

if nargin < 4
    visMode = 0;
end

% load default parameters if not set
pln = matRad_cfg.getDefaultProperties(pln,{'propOpt','propStf'});

isExternalTherapy = any(strcmp(pln.radiationMode,{'photons','protons','helium','carbon'}));
isPhotonTherapy = strcmp(pln.radiationMode,'photons');
isIonTherapy = any(strcmp(pln.radiationMode,{'protons','helium','carbon'}));
isBrachyTherapy = strcmp(pln.radiationMode,'brachy');

%% Check config and initialize
% add margin information
addmarginBool = matRad_cfg.defaults.propStf.addMargin;
if isfield(pln,'propStf') && isfield(pln.propStf,'addMargin')
    addmarginBool = pln.propStf.addMargin;
end

% get machine
try
    machine = matRad_loadMachine(pln);
catch
    matRad_cfg.dispError('Could not find the following machine file: %s',fileName);
end
%Check Config
if ~isfield(pln,'multScen')
    matRad_cfg.dispWarning('No scenario model specified! Using nominal Scenario model!');
    pln.multScen = matRad_NominalScenario(ct);
end

if isExternalTherapy
    if ~isfield(pln.propStf,'isoCenter')
        matRad_cfg.dispWarning('No isocenter specified! Using center-of-mass of all targets!');
        pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct);
    end

    if numel(pln.propStf.gantryAngles) ~= numel(pln.propStf.couchAngles)
        matRad_cfg.dispError('Inconsistent number of gantry and couch angles.');
    end

    if ~isnumeric(pln.propStf.bixelWidth) || pln.propStf.bixelWidth < 0 || ~isfinite(pln.propStf.bixelWidth)
        matRad_cfg.dispError('bixel width (spot distance) needs to be a real number [mm] larger than zero.');
    end
elseif isBrachyTherapy
    if ~isfield(pln.propStf,'needle') || ~isfield(pln.propStf.needle,'seedDistance') || ~isfield(pln.propStf.needle,'seedsNo')
        matRad_cfg.dispError('Needle information missing in pln.propStf.needle');
    end

    if ~isfield(pln.propStf,'template') || ~isfield(pln.propStf.template,'activeNeedles')
        matRad_cfg.dispError('Template information missing in pln.propStf.template!');
    end

    if ~isfield(pln.propStf,'templateRoot')
        matRad_cfg.dispError('TemplateRoot information missing in pln.propStf.templateRoot!');
    end

    if ~isa(pln.multScen,'matRad_NominalScenario') && ~strcmp(pln.multScen,'nomScen')
        matRad_cfg.dispError('Brachy Therapy does only work with a nominal scenario for now!');
    end
else
    matRad_cfg.dispError('Unknown radiation mode: %s',pln.radiationMode);
end

%Gather Initial information and computational margin extent for positioning
if isExternalTherapy
    if isIonTherapy
        if ~isfield(pln.propStf,'useRangeShifter')
            pln.propStf.useRangeShifter = false;
        end

        availableEnergies = [machine.data.energy];
        availablePeakPos  = [machine.data.peakPos] + [machine.data.offset];
        availableWidths   = [machine.data.initFocus];
        availableWidths   = [availableWidths.SisFWHMAtIso];
        maxPBwidth        = max(availableWidths) / 2.355;

        %Compute a margin to account for pencil beam width
        pbMargin = min(maxPBwidth,pln.propStf.bixelWidth);

        if pln.propStf.useRangeShifter
            %For now only a generic range shifter is used whose thickness is
            %determined by the minimum peak width to play with
            rangeShifterEqD = round(min(availablePeakPos)* 1.25);
            availablePeakPosRaShi = availablePeakPos - rangeShifterEqD;

            matRad_cfg.dispWarning('Use of range shifter enabled. matRad will generate a generic range shifter with WEPL %f to enable ranges below the shortest base data entry.',rangeShifterEqD);
        end

        if ~isfield(pln.propStf, 'longitudinalSpotSpacing')
            longitudinalSpotSpacing = matRad_cfg.propStf.defaultLongitudinalSpotSpacing;
        else
            longitudinalSpotSpacing = pln.propStf.longitudinalSpotSpacing;
        end

        if sum(availablePeakPos<0)>0
            matRad_cfg.dispError('at least one available peak position is negative - inconsistent machine file')
        end
        %clear machine;
    else
        pbMargin = pln.propStf.bixelWidth;
    end
else %Brachytherapy
    pbMargin = 0;
end

%% Prepare voxel geometry
% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln);

% find all target voxels from cst cell array
V = [];

isTarget = cellfun(@(voiType) isequal(voiType,'TARGET'),cst(:,3));
if ~any(isTarget)
    matRad_cfg.dispError('No target found in cst. Please designate at least one VOI as ''TARGET''!');
end

hasObjective = ~cellfun(@isempty,cst(:,6));
useTargetForBixelPlacement = isTarget & hasObjective;

if ~any(useTargetForBixelPlacement)
    matRad_cfg.dispWarning('No Objectives / Constraints assigned to targets. All targets will be considered for Bixel placement!');
    useTargetForBixelPlacement(isTarget) = true;
end

%Now add all used target voxels to the voxel list
for i=1:size(cst,1)
    if useTargetForBixelPlacement(i)
        V = [V; cst{i,4}{1}];
    end
end

% Remove double voxels
V = unique(V);
% generate voi cube for targets
voiTarget    = zeros(ct.cubeDim);
voiTarget(V) = 1;

%Margin info
if addmarginBool
    %Assumption for range uncertainty
    assumeRangeMargin = pln.multScen.maxAbsRangeShift + pln.multScen.maxRelRangeShift + pbMargin;

    % add margin -  account for voxel resolution, the maximum shift scenario and the current bixel width.
    margin.x  = max([ct.resolution.x max(abs(pln.multScen.isoShift(:,1)) + assumeRangeMargin)]);
    margin.y  = max([ct.resolution.y max(abs(pln.multScen.isoShift(:,2)) + assumeRangeMargin)]);
    margin.z  = max([ct.resolution.z max(abs(pln.multScen.isoShift(:,3)) + assumeRangeMargin)]);

    voiTarget = matRad_addMargin(voiTarget,cst,ct.resolution,margin,true);
    V         = find(voiTarget>0);
end

% throw error message if no target is found
if isempty(V)
    matRad_cfg.dispError('Could not find target.');
end

% get world coordinate system
ct = matRad_getWorldAxes(ct);

% Convert linear indices to 3D voxel coordinates
voxTargetWorldCoords = matRad_cubeIndex2worldCoords(V, ct);


% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(V) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

if isExternalTherapy
    % prepare structures necessary for particles
    SAD = machine.meta.SAD;

    % Define steering file like struct. Prellocating for speed.
    stf = struct;

    % loop over all angles
    for i = 1:length(pln.propStf.gantryAngles)

        % Correct for iso center position. Whit this correction Isocenter is
        % (0,0,0) [mm]
        isoCoords = voxTargetWorldCoords - pln.propStf.isoCenter(i,:);
        
        % Save meta information for treatment plan
        stf(i).gantryAngle   = pln.propStf.gantryAngles(i);
        stf(i).couchAngle    = pln.propStf.couchAngles(i);
        stf(i).bixelWidth    = pln.propStf.bixelWidth;
        stf(i).radiationMode = pln.radiationMode;
        stf(i).machine       = pln.machine;
        stf(i).SAD           = SAD;
        stf(i).isoCenter     = pln.propStf.isoCenter(i,:);

        % Get the (active) rotation matrix. We perform a passive/system
        % rotation with row vector coordinates, which would introduce two
        % inversions / transpositions of the matrix, thus no changes to the
        % rotation matrix are necessary
        rotMat_system_T = matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i));

        rot_coords = isoCoords*rotMat_system_T;

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
                for k = -floor(maxCtResolution/pln.propStf.bixelWidth):floor(maxCtResolution/pln.propStf.bixelWidth)
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
                x = [x; (x_min:pln.propStf.bixelWidth:x_max)'];
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
                stf(i).ray(j).beamletCornersAtIso = [rayPos(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    rayPos(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    rayPos(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                    rayPos(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]]*rotMat_vectors_T;
                stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                    [rayPos(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    rayPos(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                    rayPos(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                    rayPos(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMat_vectors_T;
            end
        end

        % loop over all rays to determine meta information for each ray
        stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
        
        % mm axes with isocenter at  (0,0,0)
        mmCubeIsoCenter =  matRad_world2cubeCoords(stf(i).isoCenter,ct);
        for j = stf(i).numOfRays:-1:1
            
            for ShiftScen = 1:pln.multScen.totNumShiftScen
                % ray tracing necessary to determine depth of the target
                [alphas,l{ShiftScen},rho{ShiftScen},d12,~] = matRad_siddonRayTracer(mmCubeIsoCenter + pln.multScen.isoShift(ShiftScen,:), ...
                    ct.resolution, ...
                    stf(i).sourcePoint, ...
                    stf(i).ray(j).targetPoint, ...
                    [ct.cube {voiTarget}]);

                %Used for generic range-shifter placement
                ctEntryPoint = alphas(1) * d12;
            end

            % find appropriate energies for particles
            if isIonTherapy

                % target hit
                rhoVOITarget = [];
                for shiftScen = 1:pln.multScen.totNumShiftScen
                    rhoVOITarget = [rhoVOITarget, rho{shiftScen}{end}];
                end

                if any(rhoVOITarget)
                    Counter = 0;

                    %Here we iterate through scenarios to check the required
                    %energies w.r.t lateral position.
                    %TODO: iterate over the linear scenario mask instead?
                    for CtScen = 1:pln.multScen.numOfCtScen
                        for ShiftScen = 1:pln.multScen.totNumShiftScen
                            for RangeShiftScen = 1:pln.multScen.totNumRangeScen
                                if pln.multScen.scenMask(CtScen,ShiftScen,RangeShiftScen)
                                    Counter = Counter+1;

                                    % compute radiological depths
                                    % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                                    radDepths = cumsum(l{ShiftScen} .* rho{ShiftScen}{CtScen});

                                    if pln.multScen.relRangeShift(RangeShiftScen) ~= 0 || pln.multScen.absRangeShift(RangeShiftScen) ~= 0
                                        radDepths = radDepths +...                                                        % original cube
                                            rho{ShiftScen}{CtScen}*pln.multScen.relRangeShift(RangeShiftScen) +... % rel range shift
                                            pln.multScen.absRangeShift(RangeShiftScen);                           % absolute range shift
                                        radDepths(radDepths < 0) = 0;
                                    end

                                    % find target entry & exit
                                    diff_voi    = [diff([rho{ShiftScen}{end}])];
                                    entryIx = find(diff_voi == 1);
                                    exitIx = find(diff_voi == -1);

                                    %We approximate the interface using the
                                    %rad depth between the last voxel before
                                    %and the first voxel after the interface
                                    %This captures the case that the first
                                    %relevant voxel is a target voxel
                                    targetEntry(Counter,1:length(entryIx)) = (radDepths(entryIx) + radDepths(entryIx+1)) ./ 2;
                                    targetExit(Counter,1:length(exitIx)) = (radDepths(exitIx) + radDepths(exitIx+1)) ./ 2;
                                end
                            end
                        end
                    end

                    targetEntry(targetEntry == 0) = NaN;
                    targetExit(targetExit == 0)   = NaN;

                    targetEntry = min(targetEntry);
                    targetExit  = max(targetExit);

                    %check that each energy appears only once in stf
                    if(numel(targetEntry)>1)
                        m = numel(targetEntry);
                        while(m>1)
                            if(targetEntry(m) < targetExit(m-1))
                                targetExit(m-1) = max(targetExit(m-1:m));
                                targetExit(m)=[];
                                targetEntry(m-1) = min(targetEntry(m-1:m));
                                targetEntry(m)=[];
                                m = numel(targetEntry)+1;
                            end
                            m=m-1;
                        end
                    end

                    if numel(targetEntry) ~= numel(targetExit)
                        matRad_cfg.dispError('Inconsistency during ray tracing. Please check correct assignment and overlap priorities of structure types OAR & TARGET.');
                    end

                    stf(i).ray(j).energy = [];
                    stf(i).ray(j).rangeShifter = [];

                    % Save energies in stf struct
                    for k = 1:numel(targetEntry)

                        %If we need lower energies than available, consider
                        %range shifter (if asked for)
                        if any(targetEntry < min(availablePeakPos)) && pln.propStf.useRangeShifter
                            %Get Energies to use with range shifter to fill up
                            %non-reachable low-range spots
                            raShiEnergies = availableEnergies(availablePeakPosRaShi >= targetEntry(k) & min(availablePeakPos) > availablePeakPosRaShi);

                            raShi.ID = 1;
                            raShi.eqThickness = rangeShifterEqD;
                            raShi.sourceRashiDistance = round(ctEntryPoint - 2*rangeShifterEqD,-1); %place a little away from entry, round to cms to reduce number of unique settings

                            stf(i).ray(j).energy = [stf(i).ray(j).energy raShiEnergies];
                            stf(i).ray(j).rangeShifter = [stf(i).ray(j).rangeShifter repmat(raShi,1,length(raShiEnergies))];
                        end

                        %Normal placement without rangeshifter
                        newEnergies = availableEnergies(availablePeakPos>=targetEntry(k)&availablePeakPos<=targetExit(k));


                        stf(i).ray(j).energy = [stf(i).ray(j).energy newEnergies];


                        raShi.ID = 0;
                        raShi.eqThickness = 0;
                        raShi.sourceRashiDistance = 0;
                        stf(i).ray(j).rangeShifter = [stf(i).ray(j).rangeShifter repmat(raShi,1,length(newEnergies))];
                    end


                    targetEntry = [];
                    targetExit = [];


                    % book keeping & calculate focus index
                    stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
                    currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:)',...
                        machine.meta.LUT_bxWidthminFWHM(2,:)',...
                        pln.propStf.bixelWidth, ...
                        machine.meta.LUT_bxWidthminFWHM(2,end));
                    focusIx  =  ones(stf(i).numOfBixelsPerRay(j),1);
                    [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                        repmat(stf(i).ray(j).energy,length([machine.data]),1))));

                    % get for each spot the focus index
                    for k = 1:stf(i).numOfBixelsPerRay(j)
                        focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                    end

                    stf(i).ray(j).focusIx = focusIx';

                    %Get machine bounds
                    numParticlesPerMU = 1e6*ones(1,stf(i).numOfBixelsPerRay(j));
                    minMU = zeros(1,stf(i).numOfBixelsPerRay(j));
                    maxMU = Inf(1,stf(i).numOfBixelsPerRay(j));
                    for k = 1:stf(i).numOfBixelsPerRay(j)
                        if isfield(machine.data(vEnergyIx(k)),'MUdata')
                            MUdata = machine.data(vEnergyIx(k)).MUdata;
                            if isfield(MUdata,'numParticlesPerMU')
                                numParticlesPerMU(k) = MUdata.numParticlesPerMU;
                            end

                            if isfield(MUdata,'minMU')
                                minMU(k) = MUdata.minMU;
                            end

                            if isfield(MUdata,'maxMU')
                                maxMU(k) = MUdata.maxMU;
                            end
                        end
                    end

                    stf(i).ray(j).numParticlesPerMU = numParticlesPerMU;
                    stf(i).ray(j).minMU = minMU;
                    stf(i).ray(j).maxMU = maxMU;

                else % target not hit
                    stf(i).ray(j)               = [];
                    stf(i).numOfBixelsPerRay(j) = [];
                end

            elseif isPhotonTherapy
                % book keeping for photons
                stf(i).ray(j).energy = machine.data.energy;
            else
                matRad_cfg.dispError('Error generating stf struct: invalid radiation modality.');
            end
        end

        if ~isfield(stf(i).ray,'energy')
            matRad_cfg.dispError('Error generating stf struct: no suitable energies found. Check if bixelwidth is too large.');
        end
        % store total number of rays for beam-i
        stf(i).numOfRays = size(stf(i).ray,2);

        % post processing for particle remove energy slices
        if isIonTherapy

            % get minimum energy per field
            minEnergy = min([stf(i).ray.energy]);
            maxEnergy = max([stf(i).ray.energy]);

            % get corresponding peak position
            minPeakPos  = machine.data(minEnergy == availableEnergies).peakPos;
            maxPeakPos  = machine.data(maxEnergy == availableEnergies).peakPos;

            % find set of energyies with adequate spacing


            stf(i).longitudinalSpotSpacing = longitudinalSpotSpacing;

            tolerance              = longitudinalSpotSpacing/10;

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
                ixRun = ixRun + 1;
            end

            for j = stf(i).numOfRays:-1:1
                for k = stf(i).numOfBixelsPerRay(j):-1:1
                    maskEnergy = stf(i).ray(j).energy(k) == availableEnergies;
                    if ~useEnergyBool(maskEnergy)
                        stf(i).ray(j).energy(k)         = [];
                        stf(i).ray(j).focusIx(k)        = [];
                        stf(i).ray(j).rangeShifter(k)   = [];
                        stf(i).numOfBixelsPerRay(j) = stf(i).numOfBixelsPerRay(j) - 1;
                    end
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

        % Show progress
        if matRad_cfg.logLevel > 2
            matRad_progress(i,length(pln.propStf.gantryAngles));
        end

        %% visualization
        if visMode > 0

            %center coordinates
            x = ct.x - stf(i).isoCenter(1);
            y = ct.y - stf(i).isoCenter(2);
            z = ct.z - stf(i).isoCenter(3);
            d = sqrt(sum(([x(1) y(2) z(3)] - [x(end) y(end) z(end)]).^2));
            limits = [-d/2 d/2 -d/2 d/2 -d/2 d/2];

            clf;
            % first subplot: visualization in bev
            hAxBEV = subplot(1,2,1);
            set(hAxBEV,'DataAspectRatioMode','manual');
            hold on;

            % plot rotated target coordinates
            plot3(rot_coords(:,1),rot_coords(:,2),rot_coords(:,3),'r.');



            % surface rendering
            if visMode == 2

                % computes surface
                patSurfCube      = 0*ct.cube{1};
                idx              = [cst{:,4}];
                idx              = unique(vertcat(idx{:}));
                patSurfCube(idx) = 1;

                [f,v] = isosurface(x,y,z,patSurfCube,.5);

                % rotate surface
                vRot = v*rotMat_system_T;

                % surface rendering
                surface = patch('Faces',f,'Vertices',vRot);
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
            view(hAxBEV,0,-90);
            xlabel(hAxBEV,'X [mm]');
            ylabel(hAxBEV,'Y [mm]');
            zlabel(hAxBEV,'Z [mm]');
            title(hAxBEV,'Beam''s eye view');
            
            axis(hAxBEV,limits);
            

            % second subplot: visualization in lps coordinate system
            hAxLPS = subplot(1,2,2);

            % Plot target coordinates whitout any rotation
            plot3(isoCoords(:,1),isoCoords(:,2),isoCoords(:,3),'r.');
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
            xlabel(hAxBEV,'X [mm]');
            ylabel(hAxBEV,'Y [mm]');
            zlabel(hAxBEV,'Z [mm]');
            title('lps coordinate system');
            axis(limits);
            %pause(1);
        end
        % save total number of bixels
        stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);

        % Show progress
        if matRad_cfg.logLevel > 2
            matRad_progress(i,length(pln.propStf.gantryAngles));
        end
    end
elseif isBrachyTherapy
    %translate to geometric coordinates and save in stf

    stf.targetVolume.Xvox = voxTargetWorldCoords(:,1); % angabe in mm
    stf.targetVolume.Yvox = voxTargetWorldCoords(:,2);
    stf.targetVolume.Zvox = voxTargetWorldCoords(:,3);

    %% meta info from pln
    stf.radiationMode = pln.radiationMode;
    stf.numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
    stf.numOfNeedles = nnz(pln.propStf.template.activeNeedles);
    stf.totalNumOfBixels = stf.numOfSeedsPerNeedle*stf.numOfNeedles; % means total number of seeds

    %% generate 2D template points
    % the template origin is set at its center. In the image coordinate system,
    % the center will be positioned at the bottom of the volume of interest.
    [row,col] = find(pln.propStf.template.activeNeedles);
    templX = col*pln.propStf.bixelWidth + pln.propStf.templateRoot(1) - (13+1)/2*pln.propStf.bixelWidth;
    templY = row*pln.propStf.bixelWidth + pln.propStf.templateRoot(2) - (13+1)/2*pln.propStf.bixelWidth;
    templZ = ones(size(col))                 + pln.propStf.templateRoot(3);

    stf.template = [templX';templY';templZ'];
    stf.templateNormal = [0,0,1];

    %% generate seed positions
    % seed positions can be generated from neeldes, template and oriantation
    % needles are assumed to go trough the template vertically

    % needle position
    d = pln.propStf.needle.seedDistance;
    seedsNo = pln.propStf.needle.seedsNo;
    needleDist(1,1,:) = d.*[0:seedsNo-1]'; % 1x1xN Array with seed positions on needle
    needleDir = needleDist.*[0;0;1];
    seedPos_coord_need_seed = needleDir + stf.template;
    seedPos_need_seed_coord = shiftdim(seedPos_coord_need_seed,1);
    % the output array has the dimentions (needleNo,seedNo,coordinates)
    X = seedPos_need_seed_coord(:,:,1);
    Y = seedPos_need_seed_coord(:,:,2);
    Z = seedPos_need_seed_coord(:,:,3);

    stf.seedPoints.x = reshape(X,1,[]);
    stf.seedPoints.y = reshape(Y,1,[]);
    stf.seedPoints.z = reshape(Z,1,[]);

    matRad_cfg.dispInfo('Processing completed: %d%%\n', 100);

    %%visualize results of visMode is nonzero
    % plot 3D seed positions
    if visMode > 0
        clf
        seedPoints = plot3(stf.seedPoints.x,stf.seedPoints.y,stf.seedPoints.z,'.','DisplayName', 'seed points','Color','black','markersize',5);
        title( '3D Visualization of seed points')
        xlabel('X (left) [mm]')
        ylabel('Y (posterior) [mm]')
        zlabel('Z (superior) [mm]')
        hold on

        % plot 3d VOI points
        TargX = stf.targetVolume.Xvox;
        TargY = stf.targetVolume.Yvox;
        TargZ = stf.targetVolume.Zvox;
        %Prostate = plot3(TargX,TargY,TargZ,'.', 'Color','b','DisplayName', 'prostate');

        % Prepare points for boundary calculation
        P = [TargX, TargY, TargZ];

        % Determine the environment
        if matRad_cfg.isOctave
            % Octave environment
            [uni, ~] = sort(unique(TargX));
            n = length(uni);
            outline = zeros(2*n, 2);

            for i = 1:n
                y_list = TargY(TargX == uni(i));
                y_max = max(y_list);
                y_min = min(y_list);
                outline(i, :) = [uni(i), y_max];
                outline(2*n-i+1, :) = [uni(i), y_min];
            end

            % Plot the points and the computed outline
            figure;
            plot(TargX, TargY, 'b+', 'DisplayName', 'VOI Points');
            hold on;
            plot(outline(:,1), outline(:,2), 'g-', 'LineWidth', 3, 'DisplayName', 'Computed Outline');

            % Calculate the area enclosed by the outline
            area = polyarea(outline(:,1), outline(:,2));
            disp(['Polygon area: ', num2str(area)]);

            hold off;
        else
            % MATLAB environment
            k = boundary(P, 1);
            trisurf(k, P(:,1), P(:,2), P(:,3), 'FaceColor', 'red', 'FaceAlpha', 0.1, 'LineStyle', 'none')
        end
    end

    % trow warning if seed points are more then twice the central
    % distange outsidethe TARGET volume or if no sed points are in the
    % target volume

    if (max(stf.seedPoints.x-pln.propStf.templateRoot(1)) >= 4*max(stf.targetVolume.Xvox-pln.propStf.templateRoot(1)) ||...
            min(stf.seedPoints.x-pln.propStf.templateRoot(1)) <= 4*min(stf.targetVolume.Xvox-pln.propStf.templateRoot(1)) ||...
            max(stf.seedPoints.y-pln.propStf.templateRoot(2)) >= 4*max(stf.targetVolume.Yvox-pln.propStf.templateRoot(2)) ||...
            min(stf.seedPoints.y-pln.propStf.templateRoot(2)) <= 4*min(stf.targetVolume.Yvox-pln.propStf.templateRoot(2)) || ...
            max(stf.seedPoints.z-pln.propStf.templateRoot(3)) >= 4*max(stf.targetVolume.Zvox-pln.propStf.templateRoot(3)) ||...
            min(stf.seedPoints.z-pln.propStf.templateRoot(3)) <= 4*min(stf.targetVolume.Zvox-pln.propStf.templateRoot(3)))
        matRad_cfg.dispWarning('Seeds far outside the target volume');
    end
    if (max(stf.targetVolume.Xvox) <= min(stf.seedPoints.x) || min(stf.targetVolume.Xvox) >= max(stf.seedPoints.x) ||...
            max(stf.targetVolume.Yvox) <= min(stf.seedPoints.y) || min(stf.targetVolume.Yvox) >= max(stf.seedPoints.y) ||...
            max(stf.targetVolume.Zvox) <= min(stf.seedPoints.z) || min(stf.targetVolume.Zvox) >= max(stf.seedPoints.z))
        matRad_cfg.dispWarning('no seed points in VOI')
    end
end


end
