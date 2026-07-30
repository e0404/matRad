classdef (Abstract) matRad_PhotonSequencerVMATAbstract < matRad_PhotonSequencerAbstract

    % matRad_PhotonSequencerVMATAbstract: Abstract base class for photon leaf
    %   sequencers that additionally support VMAT (dynamic/arc) delivery.
    %   Holds the machinery shared by any VMAT-capable sequencing algorithm -
    %   FMO/DAO-angle bookkeeping, arc spreading, VMAT aperture-info
    %   construction and the leaf-speed/dose-rate post-processing pipeline -
    %   which is independent of how an individual beam's shapes are decomposed.
    %   Concrete subclasses still implement sequence(); when this.runVMAT is
    %   false they behave exactly as a plain matRad_PhotonSequencerAbstract.

    properties (Constant)
        isVMATCapable = true % cheap tag for isprop-style capability checks
    end

    properties
        runVMAT = false                % toggle VMAT (dynamic/arc) delivery
        continuousAperture = false     % interpolate leaf positions between DAO angles
        arcFluenceSmoothing = 'gaussian' % fluence smoothing kernel applied before stratification, see availableArcFluenceSmoothing
        apertureSelection = 'doseAreaProduct' % how surplus apertures are discarded, see availableApertureSelection
        weightToMU = 1                 % bixel weight -> MU conversion factor (from dij.weightToMU)
    end

    methods

        function this = matRad_PhotonSequencerVMATAbstract(pln)
            if nargin < 1
                pln = [];
            end
            this = this@matRad_PhotonSequencerAbstract(pln);
        end

        function set.arcFluenceSmoothing(this, value)
            matRad_cfg = MatRad_Config.instance();
            if (ischar(value) || isstring(value)) && any(strcmp(value, this.availableArcFluenceSmoothing()))
                this.arcFluenceSmoothing = char(value);
            else
                matRad_cfg.dispError('Invalid arc fluence smoothing ''%s''! Valid values are: %s', ...
                                     char(value), strjoin(this.availableArcFluenceSmoothing(), ', '));
            end
        end

        function set.apertureSelection(this, value)
            matRad_cfg = MatRad_Config.instance();
            if (ischar(value) || isstring(value)) && any(strcmp(value, this.availableApertureSelection()))
                this.apertureSelection = char(value);
            else
                matRad_cfg.dispError('Invalid aperture selection ''%s''! Valid values are: %s', ...
                                     char(value), strjoin(this.availableApertureSelection(), ', '));
            end
        end

        function fluenceMx = smoothFluenceForArc(this, fluenceMx)
            % Fluence smoothing prior to stratification, selected by the
            % arcFluenceSmoothing property. Only applied for VMAT: it alters
            % the fluence that is decomposed, so the sequenced result no
            % longer reproduces the optimized fluence, which is exactly what
            % static/IMRT sequencing relies on.
            %
            % It defaults to 'gaussian' because VMAT needs one aperture per
            % DAO control point: the graded field rim the blur produces is
            % what lets the stratification yield more apertures as the level
            % count is raised. A perfectly flat fluence decomposes into a
            % single aperture at any number of levels, which no amount of
            % re-stratification can fix (see sequenceDynamic).

            matRad_cfg = MatRad_Config.instance();

            switch this.arcFluenceSmoothing
                case 'none'
                    % keep the optimized fluence as it is

                case 'gaussian'
                    % 5-tap Gaussian along the direction of leaf motion
                    sigma = 1;
                    kernel = exp(-((-2:2).^2) / (2 * sigma^2));
                    kernel = kernel / sum(kernel);

                    temp = zeros(size(fluenceMx));
                    for row = 1:size(fluenceMx, 1)
                        temp(row, :) = conv(fluenceMx(row, :), kernel, 'same');
                    end
                    fluenceMx = temp;

                otherwise
                    matRad_cfg.dispError('Invalid arc fluence smoothing ''%s''!', this.arcFluenceSmoothing);
            end
        end

        function sequence = applyArcSequencing(this, sequence, stf)
            % Spreads the shapes generated at each FMO-anchor beam across its
            % DAO-angle children and computes gantry rotation speed / MU rate.
            % arcSequencing needs sequence.deliveryConstraints to be set.

            if ~isfield(sequence, 'deliveryConstraints')
                sequence.deliveryConstraints = this.getMachineConstraints(stf);
            end
            sequence.beam = this.arcSequencing(sequence, stf);
        end

        function beam = arcSequencing(this, sequencing, stf)
            % Distributes the shapes generated at each FMO-anchor beam to its
            % DAO-angle children according to the sweep trajectory, then
            % computes per-beam gantry rotation speed, MU rate and rescaled
            % shape weights. Assumes shapes are already ordered in increasing
            % (left to right) leaf position.
            %
            % input:
            %   sequencing:  shape sequencing structure (needs .beam and
            %                .deliveryConstraints)
            %   stf:         matRad steering struct (beam geometry etc.)
            %
            % output:
            %   beam:        beam struct with shapes and weights distributed to
            %                the correct optimized gantry angles

            numOfBeams = numel(stf);

            beam = sequencing.beam;

            leafDir = 1;

            for i = 1:numOfBeams
                if stf(i).arc.isFMOBeam

                    % Spread apertures to each child angle
                    % according to the trajectory (mean leaf position). Assume that
                    % shapes are already order in increased (left to right) position
                    leafDir = -1 * leafDir;

                    childrenIndex = stf(i).arc.childrenIx;
                    if leafDir == -1
                        % reverse order of shapes
                        childrenIndex = flipud(childrenIndex);
                    end

                    count = 1;
                    numOfShapes = beam(i).numOfShapes;

                    for shape = 1:numOfShapes
                        childIndex = childrenIndex(count);
                        beam(childIndex).leafDir = leafDir;

                        if childIndex == i
                            % do not overwrite information, since we will need it for
                            % the remaining beams (DAO, not init)
                            beam(childIndex).tempNumOfShapes = 1;
                            beam(childIndex).tempShapes = beam(i).shapes(:, :, shape);
                            beam(childIndex).tempShapesWeight = beam(i).shapesWeight(shape);
                            beam(childIndex).fluence = beam(childIndex).tempShapes;
                            beam(childIndex).sum = beam(childIndex).tempShapesWeight * beam(childIndex).tempShapes;
                        else
                            % don't worry about overwriting
                            beam(childIndex).numOfShapes = 1;
                            beam(childIndex).shapes = beam(i).shapes(:, :, shape);
                            beam(childIndex).shapesWeight = beam(i).shapesWeight(shape);
                            beam(childIndex).fluence = beam(childIndex).shapes;
                            beam(childIndex).sum = beam(childIndex).shapesWeight * beam(childIndex).shapes;
                        end

                        count = count + 1;
                    end

                    matRad_cfg = MatRad_Config.instance();
                    sweepNames = {'right-to-left', 'left-to-right'};
                    matRad_cfg.dispDebug(['VMAT arc sequencing: FMO beam %d (%.4g deg) spread %d shape(s) ' ...
                                          'over DAO children %s, sweeping %s\n'], ...
                                         i, stf(i).gantryAngle, numOfShapes, ...
                                         mat2str(childrenIndex(1:numOfShapes)'), ...
                                         sweepNames{1 + (leafDir == 1)});
                else
                    % if beam isn't an FMO beam, then there is no info in the beam
                    % struct
                    continue
                end
            end

            for i = 1:numOfBeams
                % now go through and calculate gantry rotation speed, MU rate, etc.
                if stf(i).arc.isFMOBeam
                    beam(i).numOfShapes = beam(i).tempNumOfShapes;
                    beam(i).shapes = beam(i).tempShapes;
                    beam(i).shapesWeight = beam(i).tempShapesWeight;

                    beam(i).tempNumOfShapes = [];
                    beam(i).tempShapes = [];
                    beam(i).tempShapesWeight = [];

                    for j = 1:stf(i).arc.numOfSubChildren
                        % Prevents the aperture-info construction from attempting
                        % to convert shape to aperturevec for subchildren
                        beam(stf(i).arc.subChildrenIx(j)).numOfShapes = 0;
                    end
                end

                if stf(i).arc.isDAOBeam
                    beam(i).gantryRot = sequencing.deliveryConstraints.gantryRotationSpeed(2); % gantry rotation rate until next opt angle
                    % dose rate until next opt angle
                    beam(i).MURate = this.weightToMU .* beam(i).shapesWeight .* beam(i).gantryRot ./ ...
                        stf(i).arc.DAOAngleBordersDiff;
                    % Rescale weight to represent only this control point; weight will be shared
                    % with the interpolared control points in matRad_daoVec2ApertureInfo
                    beam(i).shapesWeight = beam(i).shapesWeight .* stf(i).arc.timeFactorCurrent;
                end
            end

            beam = rmfield(beam, {'tempShapes', 'tempShapesWeight', 'tempNumOfShapes'});

            matRad_cfg = MatRad_Config.instance();
            isDAO   = arrayfun(@(s) s.arc.isDAOBeam, stf);
            MURate  = [beam(isDAO).MURate];
            gantryRot = [beam(isDAO).gantryRot];
            matRad_cfg.dispInfo(['VMAT arc sequencing: one aperture placed on each of %d DAO control points, ' ...
                                 'MU rate %.4g-%.4g MU/s, gantry rotation %.4g deg/s\n'], ...
                                nnz(isDAO), min(MURate), max(MURate), max(gantryRot));
        end

        function constraints = getMachineConstraints(~, stf)
            % Loads the machine delivery constraints (gantry rotation speed,
            % leaf speed, MU rate) referenced by the given stf, falling back
            % to generic defaults if the machine file doesn't define any.

            matRad_cfg = MatRad_Config.instance();

            machineName = unique({stf.machine});
            radiationMode = unique({stf.radiationMode});
            if numel(machineName) > 1 || numel(radiationMode) > 1
                matRad_cfg.dispError('Mixed sequencing currently not supported for VMAT sequencing');
            end

            machine = matRad_loadMachine(struct('radiationMode', radiationMode{1}, 'machine', machineName{1}));
            if ~isfield(machine, 'constraints')
                constraints = struct( ...
                                     'gantryRotationSpeed', [0 6], ... %degree/s
                                     'leafSpeed', [0 60], ... %mm/s
                                     'monitorUnitRate', [1.25 10]); % MU/s
            else
                constraints = machine.constraints;
            end
        end

        function apertureInfo = buildVMATApertureInfo(this, sequence, stf)
            % Ports the VMAT branch of the former matRad_sequencing2ApertureInfo.m
            % to build the full apertureInfo struct including
            % apertureInfo.arc.*, then converts it to the DAO vector
            % representation via matRad_OptimizationProblemVMAT.

            % continuous aperture delivery needs the doseAngleDAO
            % bookkeeping, which the stf generator only creates when
            % pln.propSeq.continuousAperture was set before stf generation
            if this.continuousAperture
                firstDAOBeam = find(arrayfun(@(s) s.arc.isDAOBeam, stf), 1);
                if ~isempty(firstDAOBeam) && ~isfield(stf(firstDAOBeam).arc, 'doseAngleDAO')
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError(['The stf was not generated for continuous aperture delivery. ' ...
                                          'Set pln.propSeq.continuousAperture = true before matRad_generateStf.']);
                end
            end

            bixelWidth = stf(1).bixelWidth; % [mm]

            bixelIndOffset = 0;
            totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
            totalNumOfShapes = sum([sequence.beam.numOfShapes]);
            totalNumOfOptBixels = 0;
            totalNumOfLeafPairs = 0;
            vectorOffset = totalNumOfShapes + 1;
            weightOffset = 1; % used for bookkeeping in matRad_preconditionFactors
            bixOffset = 1;
            interpGetsTransition = false;

            apertureInfo.jacobiScale = ones(totalNumOfShapes, 1);
            apertureInfo.arc.jacobT = zeros(totalNumOfShapes, numel(sequence.beam));

            for i = 1:size(stf, 2)

                [geometry, bixelIndOffset] = this.getMLCGeometry(stf(i), this.numOfMLCLeafPairs, bixelIndOffset);
                dimZ = geometry.numOfActiveLeafPairs;
                minX = geometry.posOfCornerBixel(1);

                for m = 1:sequence.beam(i).numOfShapes

                    if isfield(sequence.beam(i), 'shapes')
                        shapeMap = sequence.beam(i).shapes(:, :, m);
                        leftLeafPos = NaN * ones(dimZ, 1);
                        rightLeafPos = NaN * ones(dimZ, 1);
                        for l = 1:dimZ
                            leftLeafPosInd  = find(shapeMap(l, :), 1, 'first');
                            rightLeafPosInd = find(shapeMap(l, :), 1, 'last');

                            if isempty(leftLeafPosInd) && isempty(rightLeafPosInd)
                                leftLeafPos(l) = (geometry.limLeft(l) + geometry.limRight(l)) / 2;
                                rightLeafPos(l) = leftLeafPos(l);
                            else
                                leftLeafPos(l) = (leftLeafPosInd - 1) * bixelWidth + minX - 1 / 2 * bixelWidth;
                                rightLeafPos(l) = (rightLeafPosInd - 1) * bixelWidth + minX + 1 / 2 * bixelWidth;

                                % can happen in some cases in SW trajectory sampling
                                if leftLeafPos(l) < geometry.limLeft(l)
                                    leftLeafPos(l) = geometry.limLeft(l);
                                end
                                if rightLeafPos(l) > geometry.limRight(l)
                                    rightLeafPos(l) = geometry.limRight(l);
                                end
                            end
                        end

                        apertureInfo.beam(i).shape(m).leftLeafPos = leftLeafPos;
                        apertureInfo.beam(i).shape(m).rightLeafPos = rightLeafPos;
                        apertureInfo.beam(i).shape(m).weight = sequence.beam(i).shapesWeight(m);
                        apertureInfo.beam(i).shape(m).shapeMap = shapeMap;
                    end

                    apertureInfo.beam(i).shape(m).MURate = sequence.beam(i).MURate;
                    apertureInfo.beam(i).shape(m).jacobiScale = 1;
                    apertureInfo.beam(i).shape(m).weightOffset = weightOffset;
                    weightOffset = weightOffset + 1;

                    if this.continuousAperture
                        apertureInfo.beam(i).shape(m).vectorOffset = [vectorOffset vectorOffset + dimZ];
                        vectorOffset = vectorOffset + dimZ * nnz(stf(i).arc.doseAngleDAO);
                    else
                        apertureInfo.beam(i).shape(m).vectorOffset = vectorOffset;
                        vectorOffset = vectorOffset + dimZ;
                    end
                end

                apertureInfo.beam(i).numOfShapes = sequence.beam(i).numOfShapes;
                apertureInfo.beam(i).numOfActiveLeafPairs = geometry.numOfActiveLeafPairs;
                apertureInfo.beam(i).leafPairPos = geometry.leafPairPos;
                apertureInfo.beam(i).isActiveLeafPair = geometry.isActiveLeafPair;
                apertureInfo.beam(i).centralLeafPair = geometry.centralLeafPair;
                apertureInfo.beam(i).limLeft = geometry.limLeft;
                apertureInfo.beam(i).limRight = geometry.limRight;
                apertureInfo.beam(i).bixelIndMap = geometry.bixelIndMap;
                apertureInfo.beam(i).posOfCornerBixel = geometry.posOfCornerBixel;
                apertureInfo.beam(i).MLCWindow = geometry.MLCWindow;
                apertureInfo.beam(i).gantryAngle = stf(i).gantryAngle;

                apertureInfo.beam(i).bixOffset = bixOffset;
                bixOffset = bixOffset + apertureInfo.beam(i).numOfActiveLeafPairs;

                apertureInfo.arc.beam(i).isDAOBeam = stf(i).arc.isDAOBeam;
                apertureInfo.arc.beam(i).isFMOBeam = stf(i).arc.isFMOBeam;
                apertureInfo.arc.beam(i).leafDir = sequence.beam(i).leafDir;
                apertureInfo.arc.beam(i).doseAngleBorders = stf(i).arc.doseAngleBorders;
                apertureInfo.arc.beam(i).doseAngleBorderCentreDiff = stf(i).arc.doseAngleBorderCentreDiff;
                apertureInfo.arc.beam(i).doseAngleBordersDiff = stf(i).arc.doseAngleBordersDiff;

                if apertureInfo.arc.beam(i).isDAOBeam

                    totalNumOfOptBixels = totalNumOfOptBixels + stf(i).totalNumOfBixels;
                    totalNumOfLeafPairs = totalNumOfLeafPairs + apertureInfo.beam(i).numOfShapes * apertureInfo.beam(i).numOfActiveLeafPairs;

                    apertureInfo.beam(i).gantryRot = sequence.beam(i).gantryRot;

                    apertureInfo.arc.beam(i).DAOAngleBorders = stf(i).arc.DAOAngleBorders;
                    apertureInfo.arc.beam(i).DAOAngleBorderCentreDiff = stf(i).arc.DAOAngleBorderCentreDiff;
                    apertureInfo.arc.beam(i).DAOAngleBordersDiff = stf(i).arc.DAOAngleBordersDiff;

                    apertureInfo.arc.beam(i).timeFactorCurrent = stf(i).arc.timeFactorCurrent;
                    apertureInfo.arc.beam(i).timeFactors = stf(i).arc.timeFactors;

                    apertureInfo.arc.beam(i).lastDAOBeamIx = stf(i).arc.lastDAOBeamIx;
                    apertureInfo.arc.beam(i).nextDAOBeamIx = stf(i).arc.nextDAOBeamIx;
                    apertureInfo.arc.beam(i).DAOBeamNumber = stf(i).arc.DAOBeamNumber;

                    if apertureInfo.arc.beam(i).isFMOBeam
                        apertureInfo.arc.beam(i).FMOAngleBorders = stf(i).arc.FMOAngleBorders;
                        apertureInfo.arc.beam(i).FMOAngleBorderCentreDiff = stf(i).arc.FMOAngleBorderCentreDiff;
                        apertureInfo.arc.beam(i).FMOAngleBordersDiff = stf(i).arc.FMOAngleBordersDiff;
                    end

                    if this.continuousAperture
                        apertureInfo.arc.beam(i).timeFactorIx = stf(i).arc.timeFactorIx;
                        apertureInfo.arc.beam(i).doseAngleDAO = stf(i).arc.doseAngleDAO;

                        apertureInfo.arc.beam(i).leafConstMask = 1;
                        interpGetsTransition = apertureInfo.arc.beam(i).timeFactors(3) ~= 0;
                    end

                    apertureInfo.arc.jacobT(stf(i).arc.DAOBeamNumber, i) = stf(i).arc.timeFactorCurrent;

                else
                    apertureInfo.arc.beam(i).weightFracFromLastDAO = stf(i).arc.weightFracFromLastDAO;
                    apertureInfo.arc.beam(i).timeFracFromLastDAO = stf(i).arc.timeFracFromLastDAO;
                    apertureInfo.arc.beam(i).timeFracFromNextDAO = stf(i).arc.timeFracFromNextDAO;
                    apertureInfo.arc.beam(i).lastDAOBeamIx = stf(i).arc.lastDAOBeamIx;
                    apertureInfo.arc.beam(i).nextDAOBeamIx = stf(i).arc.nextDAOBeamIx;

                    if this.continuousAperture
                        apertureInfo.arc.beam(i).weightFracFromLastDAOInitial = stf(i).arc.weightFracFromLastDAOInitial;
                        apertureInfo.arc.beam(i).weightFracFromLastDAOFinal = stf(i).arc.weightFracFromLastDAOFinal;
                        apertureInfo.arc.beam(i).weightFracFromNextDAOInitial = stf(i).arc.weightFracFromNextDAOInitial;
                        apertureInfo.arc.beam(i).weightFracFromNextDAOFinal = stf(i).arc.weightFracFromNextDAOFinal;
                    end

                    lastDAO = stf(i).arc.lastDAOBeamIx;
                    nextDAO = stf(i).arc.nextDAOBeamIx;
                    lastDoseAngleBordersDiff = stf(lastDAO).arc.doseAngleBordersDiff;

                    apertureInfo.arc.jacobT(stf(lastDAO).arc.DAOBeamNumber, i) = stf(lastDAO).arc.timeFactorCurrent .* ...
                        stf(i).arc.timeFracFromLastDAO .* stf(i).arc.doseAngleBordersDiff ./ lastDoseAngleBordersDiff;
                    apertureInfo.arc.jacobT(stf(nextDAO).arc.DAOBeamNumber, i) = stf(nextDAO).arc.timeFactorCurrent .* ...
                        stf(i).arc.timeFracFromNextDAO .* stf(i).arc.doseAngleBordersDiff ./ lastDoseAngleBordersDiff;

                    if interpGetsTransition
                        apertureInfo.arc.beam(i).leafConstMask = 1;
                    end
                    interpGetsTransition = false;
                end
            end

            % save global data
            apertureInfo.continuousAperture = this.continuousAperture;
            apertureInfo.runVMAT            = true;
            apertureInfo.preconditioner     = this.preconditioner;
            apertureInfo.bixelWidth         = bixelWidth;
            apertureInfo.numOfMLCLeafPairs  = this.numOfMLCLeafPairs;
            apertureInfo.totalNumOfBixels   = totalNumOfBixels;
            apertureInfo.totalNumOfShapes   = sum([apertureInfo.beam.numOfShapes]);
            apertureInfo.weightToMU         = this.weightToMU;

            % machine delivery constraints (reuse the ones already loaded
            % during arc sequencing when available)
            if isfield(sequence, 'deliveryConstraints')
                apertureInfo.deliveryConstraints = sequence.deliveryConstraints;
            else
                apertureInfo.deliveryConstraints = this.getMachineConstraints(stf);
            end

            apertureInfo.totalNumOfOptBixels = totalNumOfOptBixels;
            apertureInfo.doseTotalNumOfLeafPairs = sum([apertureInfo.beam(:).numOfActiveLeafPairs]);

            if apertureInfo.continuousAperture
                isDAOBeam = [apertureInfo.arc.beam.isDAOBeam];
                apertureInfo.totalNumOfLeafPairs = sum(reshape([apertureInfo.arc.beam(isDAOBeam).doseAngleDAO], 2, []), 1) * ...
                    [apertureInfo.beam(isDAOBeam).numOfActiveLeafPairs]';

                apertureInfo.arc.numLeafSpeedConstraint    = nnz([apertureInfo.arc.beam.leafConstMask]);
                apertureInfo.arc.numLeafSpeedConstraintDAO = nnz([apertureInfo.arc.beam(isDAOBeam).leafConstMask]);
            else
                apertureInfo.totalNumOfLeafPairs = totalNumOfLeafPairs;
            end

            % fix instances of leaf touching
            apertureInfo = matRad_leafTouching(apertureInfo);

            shapeInd = 0;
            for i = 1:numel(apertureInfo.beam)
                if apertureInfo.arc.beam(i).isDAOBeam
                    shapeInd = shapeInd + 1;
                    apertureInfo.arc.beam(i).timeInd = apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2 + shapeInd;
                end
            end

            % create vectors for optimization
            [apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] ...
                = matRad_OptimizationProblemVMAT.matRad_daoApertureInfo2Vec(apertureInfo);

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo(['VMAT apertureInfo: %d optimization variables ' ...
                                 '(%d shape weights + 2x%d leaf positions + %d gantry times)\n'], ...
                                numel(apertureInfo.apertureVector), apertureInfo.totalNumOfShapes, ...
                                apertureInfo.totalNumOfLeafPairs, apertureInfo.totalNumOfShapes);
            deliveryModes = {'step-and-shoot', 'continuous aperture'};
            matRad_cfg.dispInfo('VMAT apertureInfo: delivery mode %s, %d bixels total\n', ...
                                deliveryModes{1 + this.continuousAperture}, ...
                                apertureInfo.totalNumOfBixels);
            if this.continuousAperture
                matRad_cfg.dispInfo('VMAT apertureInfo: %d leaf speed constraints (%d at DAO control points)\n', ...
                                    apertureInfo.arc.numLeafSpeedConstraint, ...
                                    apertureInfo.arc.numLeafSpeedConstraintDAO);
            end
        end

        function apertureInfo = postProcessVMATApertureInfo(~, apertureInfo)
            % Interpolates sub-child gantry segments, then runs the
            % leaf-speed/dose-rate limiting pass (optDelivery).

            matRad_cfg = MatRad_Config.instance();

            apertureInfo = matRad_OptimizationProblemVMAT.matRad_daoVec2ApertureInfo(apertureInfo, apertureInfo.apertureVector);
            apertureInfo = matRad_calcMaxLeafSpeed(apertureInfo);
            leafSpeedBefore = apertureInfo.maxLeafSpeed;

            apertureInfo = matRad_enforceDeliveryConstraints(apertureInfo, 0);
            apertureInfo = matRad_calcMaxLeafSpeed(apertureInfo);

            leafSpeedLimit = apertureInfo.deliveryConstraints.leafSpeed(2);
            matRad_cfg.dispInfo(['VMAT delivery check: max leaf speed %.4g -> %.4g mm/s ' ...
                                 '(machine limit %.4g mm/s)\n'], ...
                                leafSpeedBefore, apertureInfo.maxLeafSpeed, leafSpeedLimit);
            if apertureInfo.maxLeafSpeed > leafSpeedLimit * (1 + 1e-6)
                matRad_cfg.dispWarning(['VMAT delivery check: max leaf speed %.4g mm/s still exceeds the ' ...
                                        'machine limit %.4g mm/s after slowing the gantry\n'], ...
                                       apertureInfo.maxLeafSpeed, leafSpeedLimit);
            end
        end

    end

    methods (Static)

        function smoothing = availableArcFluenceSmoothing()
            % Fluence smoothing kernels a VMAT sequencer can apply before
            % stratifying an FMO fluence map into apertures.
            smoothing = {'none', 'gaussian'};
        end

    end
end
