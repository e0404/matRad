classdef matRad_StfGeneratorPhotonVMAT < matRad_StfGeneratorPhotonRayBixelAbstract
    % matRad_StfGeneratorPhotonVMAT: STF generator for photon VMAT plans.
    %
    %   gantryAngles (inherited) are interpreted as arc anchor points: the
    %   first and last angle define the arc start/finish; any intermediate
    %   angles are waypoints the arc must pass through.  Two anchor points
    %   suffice for a simple arc, which maps cleanly to a DICOM arc export.
    %
    %   To define multiple arcs in the future, set arcIndex so that anchors
    %   belonging to the same arc share the same index value (default: 1).
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        name = 'Photon VMAT stf Generator'
        shortName = 'PhotonVMAT'
        possibleRadiationModes = {'photons'}
    end

    properties
        % Arc membership index for each anchor angle.
        % Scalar 1 (default) means all anchors belong to arc 1.
        % Set to e.g. [1 1 2 2] to define two separate arcs.
        arcIndex = 1

        % Maximum angular spacing between consecutive dose-calc beams [deg]
        maxGantryAngleSpacing    = 4
        % Maximum angular spacing between consecutive DAO control points [deg]
        maxDAOGantryAngleSpacing = 8
        % Maximum angular spacing between consecutive FMO control points [deg]
        maxFMOGantryAngleSpacing = 32

        continuousAperture = false
    end

    properties (Access = protected)
        % Internal computed angle arrays (populated by setupArcAngles).
        % These represent the full set of interpolated angles used during
        % STF generation and are not intended for direct user access.
        arcGantryAngles        % fine dose-calc angles
        arcCouchAngles         % couch angle for each fine dose-calc angle
        arcDAOGantryAngles     % direct aperture optimisation control points
        arcFMOGantryAngles     % fluence map optimisation control points
        arcStartAngle          % arc boundary start (first anchor per arc)
        arcFinishAngle         % arc boundary finish (last anchor per arc)

        % Saved user-specified anchor state, stored during initialize() so
        % that generateSourceGeometry() can restore it afterwards.
        savedAnchorGantryAngles
        savedAnchorCouchAngles
        savedIsoCenter
    end

    methods

        function this = matRad_StfGeneratorPhotonVMAT(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorPhotonRayBixelAbstract(pln);
            if isempty(this.radiationMode)
                this.radiationMode = 'photons';
            end
        end

        function setDefaults(this)
            this.setDefaults@matRad_StfGeneratorPhotonRayBixelAbstract();
            % Default to a full 360 deg arc defined by two anchor points.
            this.gantryAngles = [-180, 180];
            this.couchAngles  = [0, 0];
        end

    end

    methods (Access = protected)

        function initialize(this)
            % Override to expand isoCenter and swap in fine arc angles
            % before the base class runs, so that the base class sees the
            % correct beam count when it replicates the isoCenter.

            % Compute fine/DAO/FMO angles from the anchor points.
            this.setupArcAngles();

            nAnchors = numel(this.gantryAngles);    % current anchor count
            nFine    = numel(this.arcGantryAngles); % computed fine angle count

            % Save anchor state for restoration at end of generateSourceGeometry.
            this.savedAnchorGantryAngles = this.gantryAngles;
            this.savedAnchorCouchAngles  = this.couchAngles;
            this.savedIsoCenter          = this.isoCenter;

            % Expand isoCenter to [nFine x 3] before the base class sees it.
            % Accepted user inputs (mirroring IMRT conventions):
            %   [1 x 3]        - one isoCenter for the whole arc (same as IMRT single-iso)
            %   [nAnchors x 3] - one isoCenter per anchor point
            % Both are expanded here; the base class then only validates the size.
            if ~isempty(this.isoCenter)
                if size(this.isoCenter, 1) == 1
                    % Single isoCenter: replicate for all fine angles.
                    this.isoCenter = repmat(this.isoCenter, nFine, 1);
                elseif size(this.isoCenter, 1) == nAnchors
                    % One per anchor: assign each fine angle the isoCenter of
                    % its nearest anchor (by gantry angle).
                    isoFull = zeros(nFine, 3);
                    for k = 1:nFine
                        [~, ia] = min(abs(this.gantryAngles - this.arcGantryAngles(k)));
                        isoFull(k, :) = this.isoCenter(ia, :);
                    end
                    this.isoCenter = isoFull;
                end
                % If size is already [nFine x 3] or something else unexpected,
                % leave it untouched and let the base class validate/warn.
            end

            % Swap gantryAngles/couchAngles to fine grid so the base class
            % replicates/validates isoCenter against the correct beam count.
            this.lockAngleUpdate = true;
            this.gantryAngles    = this.arcGantryAngles;
            this.couchAngles     = this.arcCouchAngles;
            this.lockAngleUpdate = false;

            % Base class initialize: loads machine, validates/computes isoCenter,
            % builds patient geometry axes.  Fine angles are active here.
            this.initialize@matRad_StfGeneratorPhotonRayBixelAbstract();
        end

        function pbMargin = getPbMargin(this)
            pbMargin = this.bixelWidth;
        end

        function setupArcAngles(this)
            % Compute internal fine/DAO/FMO angle arrays from the user-
            % specified anchor points (this.gantryAngles) and arc grouping
            % (this.arcIndex).  Results are stored in protected properties.
            %
            % For each arc, the first and last anchor define the arc extent
            % (startingAngle / finishingAngle).  Intermediate anchors are
            % currently recorded but not yet used to subdivide the spacing
            % calculation (TODO: waypoint support).

            anchorGantry = this.gantryAngles;
            anchorCouch  = this.couchAngles;

            % Broadcast scalar arcIndex to a per-anchor vector
            if isscalar(this.arcIndex)
                arcIdx = this.arcIndex * ones(1, numel(anchorGantry));
            else
                arcIdx = this.arcIndex;
            end

            arcIds = unique(arcIdx, 'stable');

            allGantryAngles = [];
            allCouchAngles  = [];
            allDAOAngles    = [];
            allFMOAngles    = [];

            for a = 1:numel(arcIds)
                mask    = arcIdx == arcIds(a);
                anchors = anchorGantry(mask);
                couch   = anchorCouch(mask);

                startAngle  = anchors(1);
                finishAngle = anchors(end);
                couchVal    = couch(1);     % TODO: only uniform couch angle per arc

                angularRange = abs(finishAngle - startAngle);

                if this.continuousAperture
                    % In continuous mode the gantry rotates between dose
                    % positions; first/last beams are centred half a
                    % spacing inside the arc boundaries.
                    numGantryAngles    = ceil(angularRange / this.maxGantryAngleSpacing);
                    gantryAngleSpacing = angularRange / numGantryAngles;

                    numDAOGantryAngles = ceil((numGantryAngles - 1) * gantryAngleSpacing / this.maxDAOGantryAngleSpacing) + 1;
                    % Align numGantryAngles so DAO angles land exactly on fine angles
                    numGantryAngles    = (numDAOGantryAngles - 1) * ceil((numGantryAngles - 1) / (numDAOGantryAngles - 1)) + 1;
                    gantryAngleSpacing = angularRange / numGantryAngles;
                    DAOGantryAngleSpacing = (angularRange - gantryAngleSpacing) / (numDAOGantryAngles - 1);

                    firstGantryAngle = startAngle  + gantryAngleSpacing / 2;
                    lastGantryAngle  = finishAngle - gantryAngleSpacing / 2;
                else
                    % Step-and-shoot: first/last beams sit at the arc boundaries.
                    numDAOGantryAngles    = ceil(angularRange / this.maxDAOGantryAngleSpacing);
                    DAOGantryAngleSpacing = angularRange / numDAOGantryAngles;
                    numGantryAngles       = ceil(numDAOGantryAngles * DAOGantryAngleSpacing / this.maxGantryAngleSpacing);
                    gantryAngleSpacing    = angularRange / numGantryAngles;

                    firstGantryAngle = startAngle;
                    lastGantryAngle  = finishAngle;
                end

                % FMO spacing must be an odd integer multiple of the DAO spacing
                numApertures = floor(this.maxFMOGantryAngleSpacing / DAOGantryAngleSpacing);
                if mod(numApertures, 2) == 0
                    numApertures = numApertures - 1;
                end
                FMOGantryAngleSpacing = numApertures * DAOGantryAngleSpacing;
                firstFMOGantryAngle   = firstGantryAngle + DAOGantryAngleSpacing * floor(numApertures / 2);
                lastFMOGantryAngle    = lastGantryAngle  - DAOGantryAngleSpacing * floor(numApertures / 2);

                arcAngles = firstGantryAngle:gantryAngleSpacing:lastGantryAngle;
                daoAngles = firstGantryAngle:DAOGantryAngleSpacing:lastGantryAngle;
                fmoAngles = firstFMOGantryAngle:FMOGantryAngleSpacing:lastFMOGantryAngle;

                allGantryAngles = [allGantryAngles, arcAngles];
                allCouchAngles  = [allCouchAngles,  couchVal * ones(1, numel(arcAngles))];
                allDAOAngles    = [allDAOAngles,    daoAngles];
                allFMOAngles    = [allFMOAngles,    fmoAngles];
            end

            this.arcGantryAngles    = allGantryAngles;
            this.arcCouchAngles     = allCouchAngles;
            this.arcDAOGantryAngles = allDAOAngles;
            this.arcFMOGantryAngles = allFMOAngles;

            % Store arc extent boundaries for border calculations.
            % TODO: per-arc tracking when multi-arc is supported.
            this.arcStartAngle  = anchorGantry(arcIdx == arcIds(1));
            this.arcStartAngle  = this.arcStartAngle(1);
            this.arcFinishAngle = anchorGantry(arcIdx == arcIds(end));
            this.arcFinishAngle = this.arcFinishAngle(end);
        end

        function stf = generateSourceGeometry(this)
            % Fine arc angles and isoCenter are already expanded, and
            % gantryAngles/couchAngles already swapped to the fine grid by
            % initialize().  Call parent to build the per-beam stf entries.
            stf = this.generateSourceGeometry@matRad_StfGeneratorPhotonRayBixelAbstract();

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Apply VMAT configuration to stf...\n');

            %% Build master ray set: union of all per-beam ray positions,
            %  gap-filled to ensure a contiguous beam aperture.
            masterRayPosBEV = zeros(0, 3);
            for i = 1:numel(stf)
                rayPosBEV = reshape([stf(i).ray(:).rayPos_bev]', 3, stf(i).numOfRays)';
                masterRayPosBEV = union(masterRayPosBEV, rayPosBEV, 'rows');
            end

            x = masterRayPosBEV(:, 1);
            y = masterRayPosBEV(:, 2);
            z = masterRayPosBEV(:, 3);
            uniZ = unique(z);
            for j = 1:numel(uniZ)
                x_loc = x(z == uniZ(j));
                x_min = min(x_loc);
                x_max = max(x_loc);
                nNew  = (x_max - x_min) / this.bixelWidth + 1;
                x = [x; (x_min:this.bixelWidth:x_max)'];
                y = [y; zeros(nNew, 1)];
                z = [z; uniZ(j) * ones(nNew, 1)];
            end

            SAD = this.machine.meta.SAD;
            masterRayPosBEV      = unique([x, y, z], 'rows');
            masterTargetPointBEV = [2 * masterRayPosBEV(:, 1), SAD * ones(size(masterRayPosBEV, 1), 1), 2 * masterRayPosBEV(:, 3)];

            %% VMAT post-processing pass 1: assign propVMAT fields per beam
            matRad_cfg.dispInfo('VMAT stf beam type and geometry setup... ');
            stf = this.prepareArcs(stf, masterRayPosBEV,  masterTargetPointBEV);

            %% VMAT post-processing pass 2: derived quantities that require
            %  the complete propVMAT data from pass 1
            matRad_cfg.dispInfo('VMAT stf cleanup... ');

            stf = this.finalizeArcs(stf);

            % Restore object state to the user-specified anchor configuration.
            this.lockAngleUpdate = true;
            this.gantryAngles    = this.savedAnchorGantryAngles;
            this.couchAngles     = this.savedAnchorCouchAngles;
            this.lockAngleUpdate = false;
            this.isoCenter       = this.savedIsoCenter;
        end

        function stf = prepareArcs(this, stf, masterRayPosBEV,  masterTargetPointBEV)
            nBeams              = numel(stf);
            numDAO              = 1;
            DAODoseAngleBorders = zeros(2 * numel(this.arcDAOGantryAngles), 1);
            offset              = 1;
            timeFacIndOffset    = 1;
            SAD = this.machine.meta.SAD;

            for i = 1:nBeams

                %% Determine FMO parent beam
                [~, stf(i).propVMAT.beamParentFMOIndex] = min(abs(this.arcFMOGantryAngles - stf(i).gantryAngle));
                stf(i).propVMAT.beamParentGantryAngle   = this.arcFMOGantryAngles(stf(i).propVMAT.beamParentFMOIndex);
                stf(i).propVMAT.beamParentIndex         = find(abs([stf.gantryAngle] - stf(i).propVMAT.beamParentGantryAngle) < 1e-6);

                stf(i).propVMAT.FMOBeam = any(abs(this.arcFMOGantryAngles - stf(i).gantryAngle) < 1e-6);
                stf(i).propVMAT.DAOBeam = any(abs(this.arcDAOGantryAngles - stf(i).gantryAngle) < 1e-6);

                %% Dose angle borders: angular range attributed to this beam
                if i == 1
                    stf(i).propVMAT.doseAngleBorders = [this.arcStartAngle, (stf(2).gantryAngle + stf(i).gantryAngle) / 2];
                elseif i == nBeams
                    stf(i).propVMAT.doseAngleBorders = [(stf(i - 1).gantryAngle + stf(i).gantryAngle) / 2, this.arcFinishAngle];
                else
                    stf(i).propVMAT.doseAngleBorders = ([stf(i - 1).gantryAngle, stf(i + 1).gantryAngle] + stf(i).gantryAngle) / 2;
                end

                stf(i).propVMAT.doseAngleBorderCentreDiff = [stf(i).gantryAngle - stf(i).propVMAT.doseAngleBorders(1), ...
                                                             stf(i).propVMAT.doseAngleBorders(2) - stf(i).gantryAngle];
                stf(i).propVMAT.doseAngleBordersDiff = sum(stf(i).propVMAT.doseAngleBorderCentreDiff);

                if stf(i).propVMAT.DAOBeam
                    %% DAO beam: record dose angle borders and compute DAO influence range
                    DAODoseAngleBorders(offset:offset + 1) = stf(i).propVMAT.doseAngleBorders;
                    offset = offset + 2;

                    % Register as child of its FMO parent
                    parent = stf(i).propVMAT.beamParentIndex;
                    if ~isfield(stf(parent).propVMAT, 'beamChildrenGantryAngles') || isempty(stf(parent).propVMAT.beamChildrenGantryAngles)
                        stf(parent).propVMAT.numOfBeamChildren         = 0;
                        stf(parent).propVMAT.beamChildrenGantryAngles  = nan(1000, 1);
                        stf(parent).propVMAT.beamChildrenIndex         = nan(1000, 1);
                    end
                    n = stf(parent).propVMAT.numOfBeamChildren + 1;
                    stf(parent).propVMAT.numOfBeamChildren             = n;
                    stf(parent).propVMAT.beamChildrenGantryAngles(n)   = stf(i).gantryAngle;
                    stf(parent).propVMAT.beamChildrenIndex(n)          = i;

                    % DAO influence angle borders
                    DAOIndex = find(abs(this.arcDAOGantryAngles - stf(i).gantryAngle) < 1e-8);

                    if DAOIndex == 1
                        stf(i).propVMAT.DAOAngleBorders = [this.arcStartAngle, ...
                                                           (this.arcDAOGantryAngles(DAOIndex + 1) + this.arcDAOGantryAngles(DAOIndex)) / 2];
                        lastDAOIndex = i;
                        nextDAOIndex = find(abs([stf.gantryAngle] - this.arcDAOGantryAngles(DAOIndex + 1)) < 1e-8);

                    elseif DAOIndex == numel(this.arcDAOGantryAngles)
                        stf(i).propVMAT.DAOAngleBorders = [ ...
                                                           (this.arcDAOGantryAngles(DAOIndex - 1) + this.arcDAOGantryAngles(DAOIndex)) / 2, ...
                                                           this.arcFinishAngle];
                        lastDAOIndex = find(abs([stf.gantryAngle] - this.arcDAOGantryAngles(DAOIndex - 1)) < 1e-8);
                        nextDAOIndex = i;

                    else
                        stf(i).propVMAT.DAOAngleBorders = ...
                            ([this.arcDAOGantryAngles(DAOIndex - 1), this.arcDAOGantryAngles(DAOIndex + 1)] + this.arcDAOGantryAngles(DAOIndex)) / 2;
                        lastDAOIndex = i;
                        nextDAOIndex = find(abs([stf.gantryAngle] - this.arcDAOGantryAngles(DAOIndex + 1)) < 1e-8);
                    end

                    stf(i).propVMAT.lastDAOIndex = lastDAOIndex;
                    stf(i).propVMAT.nextDAOIndex = nextDAOIndex;
                    stf(i).propVMAT.DAOIndex     = numDAO;
                    numDAO = numDAO + 1;

                    stf(i).propVMAT.DAOAngleBorderCentreDiff = [stf(i).gantryAngle - stf(i).propVMAT.DAOAngleBorders(1), ...
                                                                stf(i).propVMAT.DAOAngleBorders(2) - stf(i).gantryAngle];
                    stf(i).propVMAT.DAOAngleBordersDiff = sum(stf(i).propVMAT.DAOAngleBorderCentreDiff);

                    % Time factor: fraction of DAO sector time covered by this dose sector
                    stf(i).propVMAT.timeFacCurr = stf(i).propVMAT.doseAngleBordersDiff / stf(i).propVMAT.DAOAngleBordersDiff;

                    if this.continuousAperture
                        stf(i).propVMAT.timeFac    = zeros(1, 3);
                        stf(i).propVMAT.timeFac(1) = ...
                            (stf(i).propVMAT.DAOAngleBorderCentreDiff(1) - stf(i).propVMAT.doseAngleBorderCentreDiff(1)) / ...
                            stf(i).propVMAT.DAOAngleBordersDiff;
                        stf(i).propVMAT.timeFac(2) = stf(i).propVMAT.timeFacCurr;
                        stf(i).propVMAT.timeFac(3) = ...
                            (stf(i).propVMAT.DAOAngleBorderCentreDiff(2) - stf(i).propVMAT.doseAngleBorderCentreDiff(2)) / ...
                            stf(i).propVMAT.DAOAngleBordersDiff;

                        delInd                         = stf(i).propVMAT.timeFac == 0;
                        stf(i).propVMAT.timeFacInd     = [timeFacIndOffset - 1, timeFacIndOffset, timeFacIndOffset + 1];
                        stf(i).propVMAT.timeFacInd(delInd) = 0;

                        if delInd(3)
                            timeFacIndOffset = timeFacIndOffset + 1;
                        else
                            timeFacIndOffset = timeFacIndOffset + 2;
                        end
                    else
                        stf(i).propVMAT.timeFac    = zeros(1, 2);
                        stf(i).propVMAT.timeFac(1) = stf(i).propVMAT.DAOAngleBorderCentreDiff(1) / stf(i).propVMAT.DAOAngleBordersDiff;
                        stf(i).propVMAT.timeFac(2) = stf(i).propVMAT.DAOAngleBorderCentreDiff(2) / stf(i).propVMAT.DAOAngleBordersDiff;
                    end

                else
                    %% Non-DAO beam: register as sub-child of FMO parent and record interpolation fraction
                    parent = stf(i).propVMAT.beamParentIndex;
                    if ~isfield(stf(parent).propVMAT, 'beamSubChildrenGantryAngles') || isempty(stf(parent).propVMAT.beamSubChildrenGantryAngles)
                        stf(parent).propVMAT.numOfBeamSubChildren        = 0;
                        stf(parent).propVMAT.beamSubChildrenGantryAngles = nan(1000, 1);
                        stf(parent).propVMAT.beamSubChildrenIndex        = nan(1000, 1);
                    end
                    n = stf(parent).propVMAT.numOfBeamSubChildren + 1;
                    stf(parent).propVMAT.numOfBeamSubChildren            = n;
                    stf(parent).propVMAT.beamSubChildrenGantryAngles(n)  = stf(i).gantryAngle;
                    stf(parent).propVMAT.beamSubChildrenIndex(n)         = i;

                    stf(i).propVMAT.fracFromLastDAO = (stf(nextDAOIndex).gantryAngle - stf(i).gantryAngle) / ...
                                                      (stf(nextDAOIndex).gantryAngle - stf(lastDAOIndex).gantryAngle);
                    stf(i).propVMAT.lastDAOIndex = lastDAOIndex;
                    stf(i).propVMAT.nextDAOIndex = nextDAOIndex;
                end

                %% FMO beam: compute FMO influence angle borders
                if stf(i).propVMAT.FMOBeam
                    FMOIndex = find(abs(this.arcFMOGantryAngles - stf(i).gantryAngle) < 1e-8);

                    if FMOIndex == 1
                        stf(i).propVMAT.FMOAngleBorders = [this.arcStartAngle, ...
                                                           (this.arcFMOGantryAngles(FMOIndex + 1) + this.arcFMOGantryAngles(FMOIndex)) / 2];
                    elseif FMOIndex == numel(this.arcFMOGantryAngles)
                        stf(i).propVMAT.FMOAngleBorders = [ ...
                                                           (this.arcFMOGantryAngles(FMOIndex - 1) + this.arcFMOGantryAngles(FMOIndex)) / 2, ...
                                                           this.arcFinishAngle];
                    else
                        stf(i).propVMAT.FMOAngleBorders = ...
                            ([this.arcFMOGantryAngles(FMOIndex - 1), this.arcFMOGantryAngles(FMOIndex + 1)] + this.arcFMOGantryAngles(FMOIndex)) / 2;
                    end
                    stf(i).propVMAT.FMOAngleBorderCentreDiff = [stf(i).gantryAngle - stf(i).propVMAT.FMOAngleBorders(1), ...
                                                                stf(i).propVMAT.FMOAngleBorders(2) - stf(i).gantryAngle];
                    stf(i).propVMAT.FMOAngleBordersDiff = sum(stf(i).propVMAT.FMOAngleBorderCentreDiff);
                end

                %% Assign union ray set to this beam and apply rotation
                stf(i).numOfRays         = size(masterRayPosBEV, 1);
                stf(i).numOfBixelsPerRay = ones(1, stf(i).numOfRays);
                stf(i).totalNumOfBixels  = stf(i).numOfRays;

                stf(i).sourcePoint_bev = [0, -SAD, 0];
                rotMat_vectors_T = transpose(matRad_getRotationMatrix(stf(i).gantryAngle, stf(i).couchAngle));
                stf(i).sourcePoint = stf(i).sourcePoint_bev * rotMat_vectors_T;

                for j = 1:stf(i).numOfRays
                    stf(i).ray(j).rayPos_bev      = masterRayPosBEV(j, :);
                    stf(i).ray(j).targetPoint_bev = masterTargetPointBEV(j, :);
                    stf(i).ray(j).rayPos          = masterRayPosBEV(j, :) * rotMat_vectors_T;
                    stf(i).ray(j).targetPoint     = masterTargetPointBEV(j, :) * rotMat_vectors_T;
                    stf(i).ray(j).rayCorners_SCD  = (repmat([0, this.machine.meta.SCD - SAD, 0], 4, 1) + (this.machine.meta.SCD / SAD) * ...
                                                     [masterRayPosBEV(j, :) + [+stf(i).bixelWidth / 2, 0, +stf(i).bixelWidth / 2]; ...
                                                      masterRayPosBEV(j, :) + [-stf(i).bixelWidth / 2, 0, +stf(i).bixelWidth / 2]; ...
                                                      masterRayPosBEV(j, :) + [-stf(i).bixelWidth / 2, 0, -stf(i).bixelWidth / 2]; ...
                                                      masterRayPosBEV(j, :) + [+stf(i).bixelWidth / 2, 0, -stf(i).bixelWidth / 2]]) * ...
                                                      rotMat_vectors_T;
                    stf(i).ray(j).energy = this.machine.data.energy;
                end

                matRad_progress(i, nBeams);
            end
        end

        function stf = finalizeArcs(this, stf)
            nBeams = numel(stf);
            for i = 1:nBeams
                % Remove NaN padding from child/sub-child angle lists
                if stf(i).propVMAT.FMOBeam
                    if isfield(stf(i).propVMAT, 'beamChildrenGantryAngles')
                        stf(i).propVMAT.beamChildrenGantryAngles(isnan(stf(i).propVMAT.beamChildrenGantryAngles)) = [];
                        stf(i).propVMAT.beamChildrenIndex(isnan(stf(i).propVMAT.beamChildrenIndex)) = [];
                    else
                        stf(i).propVMAT.numOfBeamChildren = 0;
                    end
                    if isfield(stf(i).propVMAT, 'beamSubChildrenGantryAngles')
                        stf(i).propVMAT.beamSubChildrenGantryAngles(isnan(stf(i).propVMAT.beamSubChildrenGantryAngles)) = [];
                        stf(i).propVMAT.beamSubChildrenIndex(isnan(stf(i).propVMAT.beamSubChildrenIndex)) = [];
                    else
                        stf(i).propVMAT.numOfBeamSubChildren = 0;
                    end
                end

                if stf(i).propVMAT.DAOBeam && this.continuousAperture
                    stf(i).propVMAT.doseAngleDAO = ones(1, 2);
                    if sum(DAODoseAngleBorders == stf(i).propVMAT.doseAngleBorders(2)) > 1
                        % Final dose angle is shared - count it only once
                        stf(i).propVMAT.doseAngleDAO(2) = 0;
                    end
                end

                if ~stf(i).propVMAT.FMOBeam && ~stf(i).propVMAT.DAOBeam
                    % Leaf position interpolation fractions
                    lastBorder = stf(stf(i).propVMAT.lastDAOIndex).propVMAT.doseAngleBorders(2);
                    nextBorder = stf(stf(i).propVMAT.nextDAOIndex).propVMAT.doseAngleBorders(1);
                    span = nextBorder - lastBorder;

                    stf(i).propVMAT.fracFromLastDAO_I = (nextBorder - stf(i).propVMAT.doseAngleBorders(1)) / span;
                    stf(i).propVMAT.fracFromLastDAO_F = (nextBorder - stf(i).propVMAT.doseAngleBorders(2)) / span;
                    stf(i).propVMAT.fracFromNextDAO_I = (stf(i).propVMAT.doseAngleBorders(1) - lastBorder) / span;
                    stf(i).propVMAT.fracFromNextDAO_F = (stf(i).propVMAT.doseAngleBorders(2) - lastBorder) / span;

                    % Time interpolation fractions (clamped to [0, 1])
                    lastDAOBorder2 = stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2);
                    stf(i).propVMAT.timeFracFromLastDAO = min(max((lastDAOBorder2 - stf(i).propVMAT.doseAngleBorders(1)) / ...
                                                                  stf(i).propVMAT.doseAngleBordersDiff, 0), 1);
                    stf(i).propVMAT.timeFracFromNextDAO = min(max((stf(i).propVMAT.doseAngleBorders(2) - lastDAOBorder2)  / ...
                                                                  stf(i).propVMAT.doseAngleBordersDiff, 0), 1);
                end

                matRad_progress(i, nBeams);
            end
        end

    end

    methods (Static)

        function [available, msg] = isAvailable(pln, machine)
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            [available, msg] = matRad_StfGeneratorPhotonRayBixelAbstract.isAvailable(pln, machine);
            if ~available
                return
            else
                available = false;
                msg = [];
            end

            try
                checkBasic    = isfield(machine, 'meta') && isfield(machine, 'data');
                checkModality = any(strcmp(matRad_StfGeneratorPhotonVMAT.possibleRadiationModes, machine.meta.radiationMode)) && ...
                                any(strcmp(matRad_StfGeneratorPhotonVMAT.possibleRadiationModes, pln.radiationMode));
                if checkModality
                    checkModality = strcmp(machine.meta.radiationMode, pln.radiationMode);
                end
                preCheck = checkBasic && checkModality;
                if ~preCheck
                    return
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic fields (meta/data/radiationMode)!';
                return
            end

            available = preCheck;
        end

    end
end
