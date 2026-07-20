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
        runVMAT = false             % toggle VMAT (dynamic/arc) delivery
        continuousAperture = false  % interpolate leaf positions between DAO angles
        weightToMU = 1              % bixel weight -> MU conversion factor (from dij.weightToMU)
    end

    methods

        function this = matRad_PhotonSequencerVMATAbstract(pln)
            if nargin < 1
                pln = [];
            end
            this = this@matRad_PhotonSequencerAbstract(pln);
        end

        function fluenceMx = smoothFluenceForArc(~, fluenceMx)
            % Gaussian fluence smoothing prior to stratification. Only used
            % for VMAT: applying this to static/IMRT sequencing would alter
            % the fluence so the sequenced result no longer reproduces the
            % optimized fluence.

            sigma = 1;
            kernel = exp(-((-2:2).^2) / (2 * sigma^2));
            kernel = kernel / sum(kernel);

            temp = zeros(size(fluenceMx));
            for row = 1:size(fluenceMx, 1)
                temp(row, :) = conv(fluenceMx(row, :), kernel, 'same');
            end
            fluenceMx = temp;
        end

        function sequence = applyArcSequencing(this, sequence, stf)
            % Spreads the shapes generated at each FMO-anchor beam across its
            % DAO-angle children and computes gantry rotation speed / MU rate.
            % arcSequencing needs sequence.constraints to be set.

            if ~isfield(sequence, 'constraints')
                sequence.constraints = this.getMachineConstraints(stf);
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
            %                .constraints)
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
                    beam(i).gantryRot = sequencing.constraints.gantryRotationSpeed(2); % gantry rotation rate until next opt angle
                    % dose rate until next opt angle
                    beam(i).MURate = this.weightToMU .* beam(i).shapesWeight .* beam(i).gantryRot ./ ...
                        stf(i).arc.DAOAngleBordersDiff;
                    % Rescale weight to represent only this control point; weight will be shared
                    % with the interpolared control points in matRad_daoVec2ApertureInfo
                    beam(i).shapesWeight = beam(i).shapesWeight .* stf(i).arc.timeFactorCurrent;
                end
            end

            beam = rmfield(beam, {'tempShapes', 'tempShapesWeight', 'tempNumOfShapes'});
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

            machine = load([radiationMode{1} '_' machineName{1}]);
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

            % machine delivery constraints
            apertureInfo.deliveryConstraints = this.getMachineConstraints(stf);

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
            apertureInfo = matRad_OptimizationProblemVMAT.leafTouching(apertureInfo);

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
        end

        function apertureInfo = postProcessVMATApertureInfo(~, apertureInfo)
            % Interpolates sub-child gantry segments, then runs the
            % leaf-speed/dose-rate limiting pass (optDelivery).

            apertureInfo = matRad_OptimizationProblemVMAT.matRad_daoVec2ApertureInfo(apertureInfo, apertureInfo.apertureVector);
            apertureInfo = matRad_OptimizationProblemVMAT.maxLeafSpeed(apertureInfo);
            apertureInfo = matRad_OptimizationProblemVMAT.optDelivery(apertureInfo, 0);
            apertureInfo = matRad_OptimizationProblemVMAT.maxLeafSpeed(apertureInfo);
        end

    end
end
