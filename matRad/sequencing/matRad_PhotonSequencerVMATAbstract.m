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
            % matRad_arcSequencing needs sequence.constraints to be set.

            if ~isfield(sequence, 'constraints')
                sequence.constraints = this.getMachineConstraints(stf);
            end
            sequence.beam = matRad_arcSequencing(sequence, stf, this.weightToMU);
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
            % apertureInfo.propVMAT.*, then converts it to the DAO vector
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
            apertureInfo.propVMAT.jacobT = zeros(totalNumOfShapes, numel(sequence.beam));

            for i = 1:size(stf, 2)

                [geometry, bixelIndOffset] = matRad_getMLCGeometry(stf(i), this.numOfMLCLeafPairs, bixelIndOffset);
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
                                leftLeafPos(l) = (geometry.lim_l(l) + geometry.lim_r(l)) / 2;
                                rightLeafPos(l) = leftLeafPos(l);
                            else
                                leftLeafPos(l) = (leftLeafPosInd - 1) * bixelWidth + minX - 1 / 2 * bixelWidth;
                                rightLeafPos(l) = (rightLeafPosInd - 1) * bixelWidth + minX + 1 / 2 * bixelWidth;

                                % can happen in some cases in SW trajectory sampling
                                if leftLeafPos(l) < geometry.lim_l(l)
                                    leftLeafPos(l) = geometry.lim_l(l);
                                end
                                if rightLeafPos(l) > geometry.lim_r(l)
                                    rightLeafPos(l) = geometry.lim_r(l);
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
                        vectorOffset = vectorOffset + dimZ * nnz(stf(i).propVMAT.doseAngleDAO);
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
                apertureInfo.beam(i).lim_l = geometry.lim_l;
                apertureInfo.beam(i).lim_r = geometry.lim_r;
                apertureInfo.beam(i).bixelIndMap = geometry.bixelIndMap;
                apertureInfo.beam(i).posOfCornerBixel = geometry.posOfCornerBixel;
                apertureInfo.beam(i).MLCWindow = geometry.MLCWindow;
                apertureInfo.beam(i).gantryAngle = stf(i).gantryAngle;

                apertureInfo.beam(i).bixOffset = bixOffset;
                bixOffset = bixOffset + apertureInfo.beam(i).numOfActiveLeafPairs;

                apertureInfo.propVMAT.beam(i).DAOBeam = stf(i).propVMAT.DAOBeam;
                apertureInfo.propVMAT.beam(i).FMOBeam = stf(i).propVMAT.FMOBeam;
                apertureInfo.propVMAT.beam(i).leafDir = sequence.beam(i).leafDir;
                apertureInfo.propVMAT.beam(i).doseAngleBorders = stf(i).propVMAT.doseAngleBorders;
                apertureInfo.propVMAT.beam(i).doseAngleBorderCentreDiff = stf(i).propVMAT.doseAngleBorderCentreDiff;
                apertureInfo.propVMAT.beam(i).doseAngleBordersDiff = stf(i).propVMAT.doseAngleBordersDiff;

                if apertureInfo.propVMAT.beam(i).DAOBeam

                    totalNumOfOptBixels = totalNumOfOptBixels + stf(i).totalNumOfBixels;
                    totalNumOfLeafPairs = totalNumOfLeafPairs + apertureInfo.beam(i).numOfShapes * apertureInfo.beam(i).numOfActiveLeafPairs;

                    apertureInfo.beam(i).gantryRot = sequence.beam(i).gantryRot;

                    apertureInfo.propVMAT.beam(i).DAOAngleBorders = stf(i).propVMAT.DAOAngleBorders;
                    apertureInfo.propVMAT.beam(i).DAOAngleBorderCentreDiff = stf(i).propVMAT.DAOAngleBorderCentreDiff;
                    apertureInfo.propVMAT.beam(i).DAOAngleBordersDiff = stf(i).propVMAT.DAOAngleBordersDiff;

                    apertureInfo.propVMAT.beam(i).timeFacCurr = stf(i).propVMAT.timeFacCurr;
                    apertureInfo.propVMAT.beam(i).timeFac = stf(i).propVMAT.timeFac;

                    apertureInfo.propVMAT.beam(i).lastDAOIndex = stf(i).propVMAT.lastDAOIndex;
                    apertureInfo.propVMAT.beam(i).nextDAOIndex = stf(i).propVMAT.nextDAOIndex;
                    apertureInfo.propVMAT.beam(i).DAOIndex = stf(i).propVMAT.DAOIndex;

                    if apertureInfo.propVMAT.beam(i).FMOBeam
                        apertureInfo.propVMAT.beam(i).FMOAngleBorders = stf(i).propVMAT.FMOAngleBorders;
                        apertureInfo.propVMAT.beam(i).FMOAngleBorderCentreDiff = stf(i).propVMAT.FMOAngleBorderCentreDiff;
                        apertureInfo.propVMAT.beam(i).FMOAngleBordersDiff = stf(i).propVMAT.FMOAngleBordersDiff;
                    end

                    if this.continuousAperture
                        apertureInfo.propVMAT.beam(i).timeFacInd = stf(i).propVMAT.timeFacInd;
                        apertureInfo.propVMAT.beam(i).doseAngleDAO = stf(i).propVMAT.doseAngleDAO;

                        apertureInfo.propVMAT.beam(i).leafConstMask = 1;
                        interpGetsTransition = apertureInfo.propVMAT.beam(i).timeFac(3) ~= 0;
                    end

                    apertureInfo.propVMAT.jacobT(stf(i).propVMAT.DAOIndex, i) = stf(i).propVMAT.timeFacCurr;

                else
                    apertureInfo.propVMAT.beam(i).fracFromLastDAO = stf(i).propVMAT.fracFromLastDAO;
                    apertureInfo.propVMAT.beam(i).timeFracFromLastDAO = stf(i).propVMAT.timeFracFromLastDAO;
                    apertureInfo.propVMAT.beam(i).timeFracFromNextDAO = stf(i).propVMAT.timeFracFromNextDAO;
                    apertureInfo.propVMAT.beam(i).lastDAOIndex = stf(i).propVMAT.lastDAOIndex;
                    apertureInfo.propVMAT.beam(i).nextDAOIndex = stf(i).propVMAT.nextDAOIndex;

                    if this.continuousAperture
                        apertureInfo.propVMAT.beam(i).fracFromLastDAO_I = stf(i).propVMAT.fracFromLastDAO_I;
                        apertureInfo.propVMAT.beam(i).fracFromLastDAO_F = stf(i).propVMAT.fracFromLastDAO_F;
                        apertureInfo.propVMAT.beam(i).fracFromNextDAO_I = stf(i).propVMAT.fracFromNextDAO_I;
                        apertureInfo.propVMAT.beam(i).fracFromNextDAO_F = stf(i).propVMAT.fracFromNextDAO_F;
                    end

                    lastDAO = stf(i).propVMAT.lastDAOIndex;
                    nextDAO = stf(i).propVMAT.nextDAOIndex;
                    lastDoseAngleBordersDiff = stf(lastDAO).propVMAT.doseAngleBordersDiff;

                    apertureInfo.propVMAT.jacobT(stf(lastDAO).propVMAT.DAOIndex, i) = stf(lastDAO).propVMAT.timeFacCurr .* ...
                        stf(i).propVMAT.timeFracFromLastDAO .* stf(i).propVMAT.doseAngleBordersDiff ./ lastDoseAngleBordersDiff;
                    apertureInfo.propVMAT.jacobT(stf(nextDAO).propVMAT.DAOIndex, i) = stf(nextDAO).propVMAT.timeFacCurr .* ...
                        stf(i).propVMAT.timeFracFromNextDAO .* stf(i).propVMAT.doseAngleBordersDiff ./ lastDoseAngleBordersDiff;

                    if interpGetsTransition
                        apertureInfo.propVMAT.beam(i).leafConstMask = 1;
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
            apertureInfo.propVMAT.constraints = this.getMachineConstraints(stf);

            apertureInfo.totalNumOfOptBixels = totalNumOfOptBixels;
            apertureInfo.doseTotalNumOfLeafPairs = sum([apertureInfo.beam(:).numOfActiveLeafPairs]);

            if apertureInfo.continuousAperture
                isDAOBeam = [apertureInfo.propVMAT.beam.DAOBeam];
                apertureInfo.totalNumOfLeafPairs = sum(reshape([apertureInfo.propVMAT.beam(isDAOBeam).doseAngleDAO], 2, []), 1) * ...
                    [apertureInfo.beam(isDAOBeam).numOfActiveLeafPairs]';

                apertureInfo.propVMAT.numLeafSpeedConstraint    = nnz([apertureInfo.propVMAT.beam.leafConstMask]);
                apertureInfo.propVMAT.numLeafSpeedConstraintDAO = nnz([apertureInfo.propVMAT.beam(isDAOBeam).leafConstMask]);
            else
                apertureInfo.totalNumOfLeafPairs = totalNumOfLeafPairs;
            end

            % fix instances of leaf touching
            apertureInfo = matRad_leafTouching(apertureInfo);

            shapeInd = 0;
            for i = 1:numel(apertureInfo.beam)
                if apertureInfo.propVMAT.beam(i).DAOBeam
                    shapeInd = shapeInd + 1;
                    apertureInfo.propVMAT.beam(i).timeInd = apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2 + shapeInd;
                end
            end

            % create vectors for optimization
            [apertureInfo.apertureVector, apertureInfo.mappingMx, apertureInfo.limMx] ...
                = matRad_OptimizationProblemVMAT.matRad_daoApertureInfo2Vec(apertureInfo);
        end

        function apertureInfo = postProcessVMATApertureInfo(~, apertureInfo)
            % Interpolates sub-child gantry segments, then runs the
            % leaf-speed/dose-rate limiting pass (matRad_optDelivery).

            apertureInfo = matRad_OptimizationProblemVMAT.matRad_daoVec2ApertureInfo(apertureInfo, apertureInfo.apertureVector);
            apertureInfo = matRad_maxLeafSpeed(apertureInfo);
            apertureInfo = matRad_optDelivery(apertureInfo, 0);
            apertureInfo = matRad_maxLeafSpeed(apertureInfo);
        end

    end
end
