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

        % Upper bounds on the realized angular spacings [deg]. They are
        % maxima, not targets: the generator picks the coarsest nesting that
        % satisfies all three, and is free to go finer where a bound or a
        % minimum count forces it to.
        % Maximum angular spacing between consecutive dose-calc beams [deg]
        maxGantryAngleSpacing    = 4
        % Maximum angular spacing between consecutive DAO control points [deg]
        maxDAOGantryAngleSpacing = 8
        % Maximum angular spacing between consecutive FMO control points [deg]
        maxFMOGantryAngleSpacing = 32

        % Minimum number of DAO control points per FMO fluence map. This is
        % the sequencer's aperture budget: each DAO control point receives
        % exactly one aperture, so a value of 1 collapses every optimized
        % fluence map to a single unmodulated field shape. Rounded up to the
        % next odd number, since the FMO angle sits at the centre of its DAO
        % children.
        minAperturesPerFMOBeam = 3

        % Optional explicit integer subdivision factors. Leave empty to
        % derive them from the max*GantryAngleSpacing bounds above; set them
        % to pin the nesting directly. Both must be odd.
        aperturesPerFMOBeam  = []   % DAO control points per FMO sector
        doseAnglesPerDAOBeam = []   % dose-calc angles per DAO sector

        % Guard against runaway problem size. Because the bounds above are
        % maxima, a tight maxFMOGantryAngleSpacing combined with a loose
        % maxDAOGantryAngleSpacing can legitimately force a very fine DAO
        % grid - fail loudly rather than silently build it.
        maxNumOfDAOAngles = 500

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

        % Guard flag set while gantryAngles/couchAngles are being
        % programmatically swapped to/from the fine angle grid.
        lockAngleUpdate = false

        % Stacked dose angle borders of all DAO beams (filled in
        % prepareArcs, consumed in finalizeArcs).
        DAODoseAngleBorders
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

        function assignPropertiesFromPln(this, pln, warnWhenPropertyChanged)
            if nargin < 3
                warnWhenPropertyChanged = false;
            end
            this.assignPropertiesFromPln@matRad_StfGeneratorPhotonRayBixelAbstract(pln, warnWhenPropertyChanged);

            % continuousAperture canonically lives under pln.propSeq (it is
            % a sequencing/delivery-mode setting that the stf generator also
            % needs); the base class only auto-maps pln.propStf.*, so bridge
            % it here. A value under pln.propStf is deprecated but honored
            % by the base-class mapping above when propSeq does not set it.
            if isfield(pln, 'propSeq') && isstruct(pln.propSeq) && isfield(pln.propSeq, 'continuousAperture')
                this.continuousAperture = pln.propSeq.continuousAperture;
            elseif isfield(pln, 'propStf') && isstruct(pln.propStf) && isfield(pln.propStf, 'continuousAperture')
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispDeprecationWarning(['pln.propStf.continuousAperture is deprecated. ' ...
                                                   'Use pln.propSeq.continuousAperture instead!']);
            end
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

            matRad_cfg = MatRad_Config.instance();
            if this.maxGantryAngleSpacing > this.maxDAOGantryAngleSpacing || ...
               this.maxDAOGantryAngleSpacing > this.maxFMOGantryAngleSpacing
                matRad_cfg.dispError(['Inconsistent VMAT gantry angle spacings: ' ...
                                      'maxGantryAngleSpacing (%g) <= maxDAOGantryAngleSpacing (%g) <= maxFMOGantryAngleSpacing (%g) required, ' ...
                                      'since dose-calc angles must subdivide DAO control points, which must subdivide FMO control points.'], ...
                                     this.maxGantryAngleSpacing, this.maxDAOGantryAngleSpacing, this.maxFMOGantryAngleSpacing);
            end

            anchorGantry = this.gantryAngles;
            anchorCouch  = this.couchAngles;

            % The photon base generator accepts a scalar couch angle and
            % applies it to every beam. VMAT expands the anchor angles
            % before the base initialize method runs, so perform the same
            % broadcasting here while the anchors are still active.
            if isscalar(anchorCouch)
                anchorCouch = repmat(anchorCouch, 1, numel(anchorGantry));
            end

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
                arcDirection = sign(finishAngle - startAngle);
                if arcDirection == 0
                    matRad_cfg.dispError('VMAT arc %g has identical start and finish angles (%g deg).', arcIds(a), startAngle);
                end

                [numFMOSectors, aperturesPerFMO, doseAnglesPerDAO] = ...
                    this.deriveArcSubdivision(angularRange, arcIds(a));

                numDAOGantryAngles = numFMOSectors * aperturesPerFMO;
                numGantryAngles    = numDAOGantryAngles * doseAnglesPerDAO;
                gantryAngleSpacing = angularRange / numGantryAngles;

                % Every beam sits at the centre of its own sector, in both
                % delivery modes. Centre sampling is what makes the three
                % angle sets nest exactly (given the odd subdivision factors)
                % and keeps every beam strictly inside the arc, so a closed
                % 360 deg arc cannot place two beams on the same physical
                % gantry position. continuousAperture changes how apertures
                % are interpolated between control points, not where the
                % control points sit.
                arcAngles = startAngle + ...
                    arcDirection * ((0:numGantryAngles - 1) + 0.5) * gantryAngleSpacing;

                % Subset by index rather than by recomputing the angles, so
                % the DAO/FMO membership tests downstream compare bit-identical
                % values instead of relying on floating-point agreement.
                daoIx = (doseAnglesPerDAO - 1) / 2 + (0:numDAOGantryAngles - 1) * doseAnglesPerDAO;
                daoAngles = arcAngles(daoIx + 1);

                fmoIx = (aperturesPerFMO - 1) / 2 + (0:numFMOSectors - 1) * aperturesPerFMO;
                fmoAngles = daoAngles(fmoIx + 1);

                allGantryAngles = [allGantryAngles, arcAngles];
                allCouchAngles  = [allCouchAngles,  couchVal * ones(1, numel(arcAngles))];
                allDAOAngles    = [allDAOAngles,    daoAngles];
                allFMOAngles    = [allFMOAngles,    fmoAngles];
            end

            this.arcGantryAngles    = allGantryAngles;
            this.arcCouchAngles     = allCouchAngles;
            this.arcDAOGantryAngles = allDAOAngles;
            this.arcFMOGantryAngles = allFMOAngles;

            if numel(arcIds) > 1
                matRad_cfg.dispInfo('VMAT total over %d arcs: %d FMO / %d DAO / %d dose angles\n', ...
                                    numel(arcIds), numel(allFMOAngles), numel(allDAOAngles), numel(allGantryAngles));
            end
            matRad_cfg.dispInfo(['VMAT: %d of %d dose angles are DAO control points, ' ...
                                 '%d are interpolated from their neighbours\n'], ...
                                numel(allDAOAngles), numel(allGantryAngles), ...
                                numel(allGantryAngles) - numel(allDAOAngles));
            matRad_cfg.dispDebug('VMAT dose angles: %s\n', mat2str(matRad_roundCompat(allGantryAngles, 4)));
            matRad_cfg.dispDebug('VMAT DAO angles:  %s\n', mat2str(matRad_roundCompat(allDAOAngles, 4)));
            matRad_cfg.dispDebug('VMAT FMO angles:  %s\n', mat2str(matRad_roundCompat(allFMOAngles, 4)));

            % Store arc extent boundaries for border calculations.
            % TODO: per-arc tracking when multi-arc is supported.
            this.arcStartAngle  = anchorGantry(arcIdx == arcIds(1));
            this.arcStartAngle  = this.arcStartAngle(1);
            this.arcFinishAngle = anchorGantry(arcIdx == arcIds(end));
            this.arcFinishAngle = this.arcFinishAngle(end);
        end

        function [numFMOSectors, aperturesPerFMO, doseAnglesPerDAO] = deriveArcSubdivision(this, angularRange, arcId)
            % Partition an arc into nested FMO / DAO / dose sectors.
            %
            % The arc holds numFMOSectors FMO sectors, each holding
            % aperturesPerFMO DAO sectors, each holding doseAnglesPerDAO dose
            % sectors. Both subdivision factors are odd so that the coarser
            % level's beam coincides with the centre one of its children.
            %
            % max*GantryAngleSpacing are upper bounds, so every count is
            % derived with ceil and every odd correction rounds *up*: going
            % finer can never violate a maximum, whereas the previous
            % floor/round-down derivation paid for the odd constraint by
            % dropping apertures, and could collapse to a single unmodulated
            % aperture per fluence map.

            matRad_cfg = MatRad_Config.instance();

            % 1. FMO sectors: coarsest partition meeting the FMO bound
            numFMOSectors = max(ceil(angularRange / this.maxFMOGantryAngleSpacing), 1);
            FMOGantryAngleSpacing = angularRange / numFMOSectors;

            % 2. DAO control points per FMO sector (= apertures per map)
            if ~isempty(this.aperturesPerFMOBeam)
                aperturesPerFMO = this.validateSubdivisionFactor(this.aperturesPerFMOBeam, 'aperturesPerFMOBeam');
            else
                aperturesPerFMO = max(ceil(FMOGantryAngleSpacing / this.maxDAOGantryAngleSpacing), ...
                                      this.minAperturesPerFMOBeam);
                aperturesPerFMO = aperturesPerFMO + mod(aperturesPerFMO + 1, 2);
            end
            DAOGantryAngleSpacing = FMOGantryAngleSpacing / aperturesPerFMO;

            % 3. dose-calc angles per DAO sector
            if ~isempty(this.doseAnglesPerDAOBeam)
                doseAnglesPerDAO = this.validateSubdivisionFactor(this.doseAnglesPerDAOBeam, 'doseAnglesPerDAOBeam');
            else
                doseAnglesPerDAO = ceil(DAOGantryAngleSpacing / this.maxGantryAngleSpacing);
                doseAnglesPerDAO = doseAnglesPerDAO + mod(doseAnglesPerDAO + 1, 2);
            end
            gantryAngleSpacing = DAOGantryAngleSpacing / doseAnglesPerDAO;

            numDAOGantryAngles = numFMOSectors * aperturesPerFMO;
            if numDAOGantryAngles > this.maxNumOfDAOAngles
                matRad_cfg.dispError(['VMAT arc %g needs %d DAO control points to satisfy the requested ' ...
                                      'spacings (max %g deg FMO / %g deg DAO) with %d apertures per fluence map, ' ...
                                      'exceeding maxNumOfDAOAngles (%d). Relax maxFMOGantryAngleSpacing, ' ...
                                      'relax maxDAOGantryAngleSpacing, or raise maxNumOfDAOAngles.'], ...
                                     arcId, numDAOGantryAngles, this.maxFMOGantryAngleSpacing, ...
                                     this.maxDAOGantryAngleSpacing, aperturesPerFMO, this.maxNumOfDAOAngles);
            end

            % Explicit subdivision factors override the bounds, so report when
            % that lets a realized spacing exceed what the user asked for.
            this.warnIfSpacingExceedsBound(gantryAngleSpacing, this.maxGantryAngleSpacing, ...
                                           'dose-calc', 'maxGantryAngleSpacing', arcId);
            this.warnIfSpacingExceedsBound(DAOGantryAngleSpacing, this.maxDAOGantryAngleSpacing, ...
                                           'DAO', 'maxDAOGantryAngleSpacing', arcId);
            this.warnIfSpacingExceedsBound(FMOGantryAngleSpacing, this.maxFMOGantryAngleSpacing, ...
                                           'FMO', 'maxFMOGantryAngleSpacing', arcId);

            if aperturesPerFMO < 3
                matRad_cfg.dispWarning(['VMAT arc %g yields only %d aperture(s) per FMO fluence map. ' ...
                                        'Each DAO control point receives exactly one aperture, so the ' ...
                                        'sequenced plan will carry little or no modulation.'], ...
                                       arcId, aperturesPerFMO);
            end

            matRad_cfg.dispInfo(['VMAT arc %g: %d FMO / %d DAO / %d dose angles, %d apertures per fluence map, ' ...
                                 'spacings %.3g / %.3g / %.3g deg (FMO/DAO/dose)\n'], ...
                                arcId, numFMOSectors, numDAOGantryAngles, ...
                                numDAOGantryAngles * doseAnglesPerDAO, aperturesPerFMO, ...
                                FMOGantryAngleSpacing, DAOGantryAngleSpacing, gantryAngleSpacing);
        end

        function warnIfSpacingExceedsBound(~, realized, bound, levelName, propName, arcId)
            if realized > bound * (1 + 1e-12)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning(['VMAT arc %g: realized %s angle spacing (%g deg) exceeds %s (%g deg) ' ...
                                        'because the subdivision factors were set explicitly.'], ...
                                       arcId, levelName, realized, propName, bound);
            end
        end

        function value = validateSubdivisionFactor(~, value, propName)
            matRad_cfg = MatRad_Config.instance();
            if ~isscalar(value) || ~isfinite(value) || value < 1 || mod(value, 1) ~= 0
                matRad_cfg.dispError('%s must be a positive integer scalar, got %s.', ...
                                     propName, mat2str(value));
            end
            if mod(value, 2) == 0
                matRad_cfg.dispError(['%s must be odd (got %d) so that the coarser control point ' ...
                                      'coincides with the centre one of its children.'], propName, value);
            end
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

            %% VMAT post-processing pass 1: assign arc fields per beam
            matRad_cfg.dispInfo('VMAT stf beam type and geometry setup... ');
            stf = this.prepareArcs(stf, masterRayPosBEV,  masterTargetPointBEV);

            %% VMAT post-processing pass 2: derived quantities that require
            %  the complete arc data from pass 1
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
            this.DAODoseAngleBorders = zeros(2 * numel(this.arcDAOGantryAngles), 1);
            offset              = 1;
            timeFacIndOffset    = 1;
            SAD = this.machine.meta.SAD;

            % Beams sit at the centre of their dose sector, so the first and
            % last dose beams of an arc lie *outside* the span bracketed by
            % DAO control points. Resolve each beam's bracketing DAO beams up
            % front rather than carrying them forward from the previous loop
            % iteration, which silently assumed beam 1 was a DAO beam.
            daoBeamIx = find(arrayfun(@(s) any(abs(this.arcDAOGantryAngles - s.gantryAngle) < 1e-6), stf));

            for i = 1:nBeams

                %% Determine FMO parent beam
                [~, stf(i).arc.parentFMOIx] = min(abs(this.arcFMOGantryAngles - stf(i).gantryAngle));
                stf(i).arc.parentGantryAngle   = this.arcFMOGantryAngles(stf(i).arc.parentFMOIx);
                [~, stf(i).arc.parentIx]    = min(abs([stf.gantryAngle] - stf(i).arc.parentGantryAngle));

                stf(i).arc.isFMOBeam = any(abs(this.arcFMOGantryAngles - stf(i).gantryAngle) < 1e-6);
                stf(i).arc.isDAOBeam = any(abs(this.arcDAOGantryAngles - stf(i).gantryAngle) < 1e-6);

                %% Dose angle borders: angular range attributed to this beam
                if i == 1
                    stf(i).arc.doseAngleBorders = [this.arcStartAngle, (stf(2).gantryAngle + stf(i).gantryAngle) / 2];
                elseif i == nBeams
                    stf(i).arc.doseAngleBorders = [(stf(i - 1).gantryAngle + stf(i).gantryAngle) / 2, this.arcFinishAngle];
                else
                    stf(i).arc.doseAngleBorders = ([stf(i - 1).gantryAngle, stf(i + 1).gantryAngle] + stf(i).gantryAngle) / 2;
                end

                stf(i).arc.doseAngleBorderCentreDiff = abs([stf(i).gantryAngle - stf(i).arc.doseAngleBorders(1), ...
                                                            stf(i).arc.doseAngleBorders(2) - stf(i).gantryAngle]);
                stf(i).arc.doseAngleBordersDiff = sum(stf(i).arc.doseAngleBorderCentreDiff);

                if stf(i).arc.isDAOBeam
                    %% DAO beam: record dose angle borders and compute DAO influence range
                    this.DAODoseAngleBorders(offset:offset + 1) = stf(i).arc.doseAngleBorders;
                    offset = offset + 2;

                    % Register as child of its FMO parent
                    parent = stf(i).arc.parentIx;
                    if ~isfield(stf(parent).arc, 'childrenGantryAngles') || isempty(stf(parent).arc.childrenGantryAngles)
                        stf(parent).arc.numOfChildren         = 0;
                        stf(parent).arc.childrenGantryAngles  = nan(1000, 1);
                        stf(parent).arc.childrenIx         = nan(1000, 1);
                    end
                    n = stf(parent).arc.numOfChildren + 1;
                    stf(parent).arc.numOfChildren             = n;
                    stf(parent).arc.childrenGantryAngles(n)   = stf(i).gantryAngle;
                    stf(parent).arc.childrenIx(n)          = i;

                    % DAO influence angle borders
                    DAOBeamNumber   = find(daoBeamIx == i);
                    numDAOBeams     = numel(daoBeamIx);
                    neighbourAngles = [NaN, NaN];
                    if DAOBeamNumber > 1
                        neighbourAngles(1) = stf(daoBeamIx(DAOBeamNumber - 1)).gantryAngle;
                    end
                    if DAOBeamNumber < numDAOBeams
                        neighbourAngles(2) = stf(daoBeamIx(DAOBeamNumber + 1)).gantryAngle;
                    end

                    % the outermost DAO sectors run out to the arc boundary
                    stf(i).arc.DAOAngleBorders = (neighbourAngles + stf(i).gantryAngle) / 2;
                    if DAOBeamNumber == 1
                        stf(i).arc.DAOAngleBorders(1) = this.arcStartAngle;
                    end
                    if DAOBeamNumber == numDAOBeams
                        stf(i).arc.DAOAngleBorders(2) = this.arcFinishAngle;
                    end

                    % a DAO beam anchors the interval starting at it, except
                    % the last one, which anchors the interval ending at it
                    if DAOBeamNumber < numDAOBeams
                        lastDAOBeamIx = i;
                        nextDAOBeamIx = daoBeamIx(DAOBeamNumber + 1);
                    else
                        lastDAOBeamIx = daoBeamIx(max(DAOBeamNumber - 1, 1));
                        nextDAOBeamIx = i;
                    end

                    stf(i).arc.lastDAOBeamIx = lastDAOBeamIx;
                    stf(i).arc.nextDAOBeamIx = nextDAOBeamIx;
                    stf(i).arc.DAOBeamNumber     = numDAO;
                    numDAO = numDAO + 1;

                    stf(i).arc.DAOAngleBorderCentreDiff = abs([stf(i).gantryAngle - stf(i).arc.DAOAngleBorders(1), ...
                                                               stf(i).arc.DAOAngleBorders(2) - stf(i).gantryAngle]);
                    stf(i).arc.DAOAngleBordersDiff = sum(stf(i).arc.DAOAngleBorderCentreDiff);

                    % Time factor: fraction of DAO sector time covered by this dose sector
                    stf(i).arc.timeFactorCurrent = stf(i).arc.doseAngleBordersDiff / stf(i).arc.DAOAngleBordersDiff;

                    if this.continuousAperture
                        stf(i).arc.timeFactors    = zeros(1, 3);
                        stf(i).arc.timeFactors(1) = ...
                            (stf(i).arc.DAOAngleBorderCentreDiff(1) - stf(i).arc.doseAngleBorderCentreDiff(1)) / ...
                            stf(i).arc.DAOAngleBordersDiff;
                        stf(i).arc.timeFactors(2) = stf(i).arc.timeFactorCurrent;
                        stf(i).arc.timeFactors(3) = ...
                            (stf(i).arc.DAOAngleBorderCentreDiff(2) - stf(i).arc.doseAngleBorderCentreDiff(2)) / ...
                            stf(i).arc.DAOAngleBordersDiff;

                        % timeFactors(1) and (3) index the leaf-sweep segments
                        % shared with the previous and next DAO control point.
                        % The outermost DAO sectors run out to the arc boundary,
                        % where there is no neighbouring aperture to sweep to -
                        % the dose beams out there simply repeat this aperture -
                        % so those segments do not exist. Zeroing them keeps
                        % timeFactorIx (which numbers the segments, and starts at
                        % timeFacIndOffset - 1) from producing a 0 index.
                        if DAOBeamNumber == 1
                            stf(i).arc.timeFactors(1) = 0;
                        end
                        if DAOBeamNumber == numDAOBeams
                            stf(i).arc.timeFactors(3) = 0;
                        end

                        delInd                         = stf(i).arc.timeFactors == 0;
                        stf(i).arc.timeFactorIx     = [timeFacIndOffset - 1, timeFacIndOffset, timeFacIndOffset + 1];
                        stf(i).arc.timeFactorIx(delInd) = 0;

                        if delInd(3)
                            timeFacIndOffset = timeFacIndOffset + 1;
                        else
                            timeFacIndOffset = timeFacIndOffset + 2;
                        end
                    else
                        stf(i).arc.timeFactors    = zeros(1, 2);
                        stf(i).arc.timeFactors(1) = stf(i).arc.DAOAngleBorderCentreDiff(1) / stf(i).arc.DAOAngleBordersDiff;
                        stf(i).arc.timeFactors(2) = stf(i).arc.DAOAngleBorderCentreDiff(2) / stf(i).arc.DAOAngleBordersDiff;
                    end

                else
                    %% Non-DAO beam: register as sub-child of FMO parent and record interpolation fraction
                    parent = stf(i).arc.parentIx;
                    if ~isfield(stf(parent).arc, 'subChildrenGantryAngles') || isempty(stf(parent).arc.subChildrenGantryAngles)
                        stf(parent).arc.numOfSubChildren        = 0;
                        stf(parent).arc.subChildrenGantryAngles = nan(1000, 1);
                        stf(parent).arc.subChildrenIx        = nan(1000, 1);
                    end
                    n = stf(parent).arc.numOfSubChildren + 1;
                    stf(parent).arc.numOfSubChildren            = n;
                    stf(parent).arc.subChildrenGantryAngles(n)  = stf(i).gantryAngle;
                    stf(parent).arc.subChildrenIx(n)         = i;

                    lastDAOBeamIx = daoBeamIx(find(daoBeamIx < i, 1, 'last'));
                    nextDAOBeamIx = daoBeamIx(find(daoBeamIx > i, 1, 'first'));
                    if isempty(lastDAOBeamIx)
                        lastDAOBeamIx = nextDAOBeamIx;
                    end
                    if isempty(nextDAOBeamIx)
                        nextDAOBeamIx = lastDAOBeamIx;
                    end

                    if lastDAOBeamIx == nextDAOBeamIx
                        % leading/trailing dose beam: it lies inside the first
                        % (last) DAO sector but before (after) that sector's
                        % control point, so its aperture comes entirely from
                        % that single control point
                        stf(i).arc.weightFracFromLastDAO = 1;
                    else
                        stf(i).arc.weightFracFromLastDAO = (stf(nextDAOBeamIx).gantryAngle - stf(i).gantryAngle) / ...
                                                          (stf(nextDAOBeamIx).gantryAngle - stf(lastDAOBeamIx).gantryAngle);
                    end
                    stf(i).arc.lastDAOBeamIx = lastDAOBeamIx;
                    stf(i).arc.nextDAOBeamIx = nextDAOBeamIx;
                end

                %% FMO beam: compute FMO influence angle borders
                if stf(i).arc.isFMOBeam
                    FMOIndex = find(abs(this.arcFMOGantryAngles - stf(i).gantryAngle) < 1e-8);

                    if FMOIndex == 1
                        stf(i).arc.FMOAngleBorders = [this.arcStartAngle, ...
                                                      (this.arcFMOGantryAngles(FMOIndex + 1) + this.arcFMOGantryAngles(FMOIndex)) / 2];
                    elseif FMOIndex == numel(this.arcFMOGantryAngles)
                        stf(i).arc.FMOAngleBorders = [ ...
                                                      (this.arcFMOGantryAngles(FMOIndex - 1) + this.arcFMOGantryAngles(FMOIndex)) / 2, ...
                                                      this.arcFinishAngle];
                    else
                        stf(i).arc.FMOAngleBorders = ...
                            ([this.arcFMOGantryAngles(FMOIndex - 1), this.arcFMOGantryAngles(FMOIndex + 1)] + this.arcFMOGantryAngles(FMOIndex)) / 2;
                    end
                    stf(i).arc.FMOAngleBorderCentreDiff = abs([stf(i).gantryAngle - stf(i).arc.FMOAngleBorders(1), ...
                                                               stf(i).arc.FMOAngleBorders(2) - stf(i).gantryAngle]);
                    stf(i).arc.FMOAngleBordersDiff = sum(stf(i).arc.FMOAngleBorderCentreDiff);
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
                if stf(i).arc.isFMOBeam
                    if isfield(stf(i).arc, 'childrenGantryAngles')
                        stf(i).arc.childrenGantryAngles(isnan(stf(i).arc.childrenGantryAngles)) = [];
                        stf(i).arc.childrenIx(isnan(stf(i).arc.childrenIx)) = [];
                    else
                        stf(i).arc.numOfChildren = 0;
                    end
                    if isfield(stf(i).arc, 'subChildrenGantryAngles')
                        stf(i).arc.subChildrenGantryAngles(isnan(stf(i).arc.subChildrenGantryAngles)) = [];
                        stf(i).arc.subChildrenIx(isnan(stf(i).arc.subChildrenIx)) = [];
                    else
                        stf(i).arc.numOfSubChildren = 0;
                    end
                end

                if stf(i).arc.isDAOBeam && this.continuousAperture
                    stf(i).arc.doseAngleDAO = ones(1, 2);
                    if sum(this.DAODoseAngleBorders == stf(i).arc.doseAngleBorders(2)) > 1
                        % Final dose angle is shared - count it only once
                        stf(i).arc.doseAngleDAO(2) = 0;
                    end
                end

                if ~stf(i).arc.isFMOBeam && ~stf(i).arc.isDAOBeam
                    if stf(i).arc.lastDAOBeamIx == stf(i).arc.nextDAOBeamIx
                        % leading/trailing dose beam - no pair to interpolate
                        % between, everything comes from the one DAO beam
                        % whose sector it sits in (see prepareArcs)
                        stf(i).arc.weightFracFromLastDAOInitial = 1;
                        stf(i).arc.weightFracFromLastDAOFinal   = 1;
                        stf(i).arc.weightFracFromNextDAOInitial = 0;
                        stf(i).arc.weightFracFromNextDAOFinal   = 0;
                        stf(i).arc.timeFracFromLastDAO          = 1;
                        stf(i).arc.timeFracFromNextDAO          = 0;
                        matRad_progress(i, nBeams);
                        continue
                    end

                    % Leaf position interpolation fractions
                    lastBorder = stf(stf(i).arc.lastDAOBeamIx).arc.doseAngleBorders(2);
                    nextBorder = stf(stf(i).arc.nextDAOBeamIx).arc.doseAngleBorders(1);
                    span = nextBorder - lastBorder;

                    stf(i).arc.weightFracFromLastDAOInitial = (nextBorder - stf(i).arc.doseAngleBorders(1)) / span;
                    stf(i).arc.weightFracFromLastDAOFinal = (nextBorder - stf(i).arc.doseAngleBorders(2)) / span;
                    stf(i).arc.weightFracFromNextDAOInitial = (stf(i).arc.doseAngleBorders(1) - lastBorder) / span;
                    stf(i).arc.weightFracFromNextDAOFinal = (stf(i).arc.doseAngleBorders(2) - lastBorder) / span;

                    % Time interpolation fractions (clamped to [0, 1]).
                    % The DAO border splits this dose sector into a part
                    % belonging to the last and one belonging to the next DAO
                    % beam, so the two fractions must sum to 1. Keep both
                    % numerator and denominator signed (as for the weight
                    % fractions above): the arc direction then cancels, while
                    % a DAO border lying outside the sector still clamps to
                    % 0 / 1 rather than to 1 / 1.
                    lastDAOBorder2 = stf(stf(i).arc.lastDAOBeamIx).arc.DAOAngleBorders(2);
                    doseSectorSpan = stf(i).arc.doseAngleBorders(2) - stf(i).arc.doseAngleBorders(1);
                    stf(i).arc.timeFracFromLastDAO = min(max((lastDAOBorder2 - stf(i).arc.doseAngleBorders(1)) / ...
                                                             doseSectorSpan, 0), 1);
                    stf(i).arc.timeFracFromNextDAO = min(max((stf(i).arc.doseAngleBorders(2) - lastDAOBorder2) / ...
                                                             doseSectorSpan, 0), 1);
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
