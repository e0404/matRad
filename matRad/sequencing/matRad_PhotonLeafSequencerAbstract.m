classdef  (Abstract) matRad_PhotonLeafSequencerAbstract < matRad_SequencerBase

    % matRad_PhotonLeafSequencerAbstract: Abstract base class for photon multileaf
    %   collimator (MLC) leaf sequencers. Provides the shared fluence
    %   stratification, aperture-info generation and segment visualization;
    %   concrete subclasses implement the specific leaf sequencing algorithm.
    properties
        numOfMLCLeafPairs = 80
        numLevels = 5 % number of fluence stratification levels for leaf sequencing (for VMAT only a
        % starting value: it is increased automatically until every FMO beam yields
        % at least as many apertures as it has DAO child angles)
        preconditioner = false % apply matRad_preconditionFactors to the resulting apertureInfo
    end

    % Deprecated properties referencing a newer one
    properties (Dependent)
        sequencingLevel % deprecated alias for numLevels
    end

    methods

        function this = matRad_PhotonLeafSequencerAbstract(pln)
            % Constructor, forwards to the base class. An explicit constructor
            % chain is required so that property assignment from pln persists
            % under Octave.
            if nargin < 1
                pln = [];
            end
            this = this@matRad_SequencerBase(pln);
        end

        function set.sequencingLevel(this, value)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The sequencer property sequencingLevel is deprecated. Use numLevels instead!');
            this.numLevels = value;
        end

        function value = get.sequencingLevel(this)
            value = this.numLevels;
        end

        function sequence = sequence(this, w, stf)

            throw(MException('MATLAB:class:AbstractMember', 'Abstract function sequence needs to be implemented!'));
        end

        function [d0, dCurrent, shapes, calFac, indInMx] = initBeam(this, stf, wCurr)

            numOfRaysPerBeam = size(stf.ray, 2);
            X = ones(numOfRaysPerBeam, 1) * NaN;
            Z = ones(numOfRaysPerBeam, 1) * NaN;

            for j = 1:numOfRaysPerBeam
                X(j) = stf.ray(j).rayPos_bev(:, 1);
                Z(j) = stf.ray(j).rayPos_bev(:, 3);
            end

            % sort bixels into matrix
            minX = min(X);
            maxX = max(X);
            minZ = min(Z);
            maxZ = max(Z);

            dimOfFluenceMxX = (maxX - minX) / stf.bixelWidth + 1;
            dimOfFluenceMxZ = (maxZ - minZ) / stf.bixelWidth + 1;

            % Create the fluence matrix.
            fluenceMx = zeros(dimOfFluenceMxZ, dimOfFluenceMxX);

            % Calculate X and Z position of every fluence's matrix spot z axis =
            % axis of leaf movement!
            xPos = (X - minX) / stf.bixelWidth + 1;
            zPos = (Z - minZ) / stf.bixelWidth + 1;

            % Make subscripts for fluence matrix
            indInMx = zPos + (xPos - 1) * dimOfFluenceMxZ;

            % Save weights in fluence matrix.
            fluenceMx(indInMx) = wCurr .* ones(numOfRaysPerBeam, 1);

            % Stratification
            calFac = max(fluenceMx(:));
            dCurrent = round(fluenceMx / calFac * this.numLevels);

            % Save the stratification in the initial intensity matrix d0.
            d0 = dCurrent;

            % container to remember generated shapes; allocate space for 10000 shapes
            shapes = NaN * ones(dimOfFluenceMxZ, dimOfFluenceMxX, 10000);

        end

        function sequence = sequencing2ApertureInfo(this, sequence, stf)
            % MLC parameters:
            bixelWidth = stf(1).bixelWidth; % [mm]

            % initializing variables
            bixelIndOffset = 0; % used for creation of bixel index maps
            totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
            totalNumOfShapes = sum([sequence.beam.numOfShapes]);
            vectorOffset = totalNumOfShapes + 1; % used for bookkeeping in the vector for optimization
            weightOffset = 1; % used for bookkeeping in matRad_preconditionFactors

            % loop over all beams
            for i = 1:size(stf, 2)

                %% 1. read stf and derive the MLC geometry (Ray & Bixelindex maps)
                [geometry, bixelIndOffset] = this.getMLCGeometry(stf(i), this.numOfMLCLeafPairs, bixelIndOffset);
                dimZ = geometry.numOfActiveLeafPairs;
                minX = geometry.posOfCornerBixel(1);

                % get leaf positions for all shapes
                % leaf positions can be extracted from the shapes created in Sequencing
                for m = 1:sequence.beam(i).numOfShapes

                    % loading shape from Sequencing result
                    shapeMap = sequence.beam(i).shapes(:, :, m);
                    % get left and right leaf indices from shapemap
                    % initializing limits
                    leftLeafPos = NaN * ones(dimZ, 1);
                    rightLeafPos = NaN * ones(dimZ, 1);
                    % looping over leaf pairs
                    for l = 1:dimZ
                        leftLeafPosInd  = find(shapeMap(l, :), 1, 'first');
                        rightLeafPosInd = find(shapeMap(l, :), 1, 'last');

                        if isempty(leftLeafPosInd) && isempty(rightLeafPosInd) % if no bixel is open, use limits from Ray positions
                            leftLeafPos(l) = (geometry.limLeft(l) + geometry.limRight(l)) / 2;
                            rightLeafPos(l) = leftLeafPos(l);
                        else
                            % the physical position [mm] can be calculated from the indices
                            leftLeafPos(l) = (leftLeafPosInd - 1) * bixelWidth + ...
                                minX - 1 / 2 * bixelWidth;
                            rightLeafPos(l) = (rightLeafPosInd - 1) * bixelWidth + ...
                                minX + 1 / 2 * bixelWidth;

                        end
                    end

                    % save data for each shape of this beam
                    sequence.apertureInfo.beam(i).shape(m).leftLeafPos = leftLeafPos;
                    sequence.apertureInfo.beam(i).shape(m).rightLeafPos = rightLeafPos;
                    sequence.apertureInfo.beam(i).shape(m).weight = sequence.beam(i).shapesWeight(m);
                    sequence.apertureInfo.beam(i).shape(m).shapeMap = shapeMap;
                    sequence.apertureInfo.beam(i).shape(m).vectorOffset = vectorOffset;
                    sequence.apertureInfo.beam(i).shape(m).jacobiScale = 1;
                    sequence.apertureInfo.beam(i).shape(m).weightOffset = weightOffset;

                    % update index for bookkeeping
                    vectorOffset = vectorOffset + dimZ;
                    weightOffset = weightOffset + 1;

                end

                % save data for each beam
                sequence.apertureInfo.beam(i).numOfShapes = sequence.beam(i).numOfShapes;
                sequence.apertureInfo.beam(i).numOfActiveLeafPairs = geometry.numOfActiveLeafPairs;
                sequence.apertureInfo.beam(i).leafPairPos = geometry.leafPairPos;
                sequence.apertureInfo.beam(i).isActiveLeafPair = geometry.isActiveLeafPair;
                sequence.apertureInfo.beam(i).centralLeafPair = geometry.centralLeafPair;
                sequence.apertureInfo.beam(i).limLeft = geometry.limLeft;
                sequence.apertureInfo.beam(i).limRight = geometry.limRight;
                sequence.apertureInfo.beam(i).bixelIndMap = geometry.bixelIndMap;
                sequence.apertureInfo.beam(i).posOfCornerBixel = geometry.posOfCornerBixel;
                sequence.apertureInfo.beam(i).MLCWindow = geometry.MLCWindow;

            end

            % save global data
            sequence.apertureInfo.runVMAT = false;
            sequence.apertureInfo.preconditioner = this.preconditioner;
            sequence.apertureInfo.bixelWidth = bixelWidth;
            sequence.apertureInfo.numOfMLCLeafPairs = this.numOfMLCLeafPairs;
            sequence.apertureInfo.totalNumOfBixels = totalNumOfBixels;
            sequence.apertureInfo.totalNumOfShapes = sum([sequence.apertureInfo.beam.numOfShapes]);
            sequence.apertureInfo.totalNumOfLeafPairs = sum([sequence.apertureInfo.beam.numOfShapes] * ...
                                                            [sequence.apertureInfo.beam.numOfActiveLeafPairs]');
            sequence.apertureInfo.doseTotalNumOfLeafPairs = sequence.apertureInfo.totalNumOfLeafPairs;
            sequence.apertureInfo.jacobiScale = ones(sequence.apertureInfo.totalNumOfShapes, 1);

            % create vectors for optimization
            [sequence.apertureInfo.apertureVector, sequence.apertureInfo.mappingMx, sequence.apertureInfo.limMx] ...
                = matRad_OptimizationProblemDAO.matRad_daoApertureInfo2Vec(sequence.apertureInfo);
        end

        function plotSegments(this, sequencing)
            % create the sequencing figure
            sz = [800 1000]; % figure size
            screensize = get(0, 'ScreenSize');
            xpos = ceil((screensize(3) - sz(2)) / 2); % center the figure on the screen horizontally
            ypos = ceil((screensize(4) - sz(1)) / 2); % center the figure on the screen vertically
            seqFig = figure('position', [xpos, ypos, sz(2), sz(1)]);

            for i = 1:numel(sequencing.beam)

                D_0 = sequencing.beam(i).fluence;

                clf(seqFig);
                colormap(seqFig, 'jet');

                seqSubPlots(1) = subplot(2, 2, 1, 'parent', seqFig);
                imagesc(sequencing.beam(i).fluence, 'parent', seqSubPlots(1));
                set(seqSubPlots(1), 'CLim', [0 this.numLevels], 'YDir', 'normal');
                title(seqSubPlots(1), ['Beam # ' num2str(i) ': max(D_0) = ' num2str(max(D_0(:))) ...
                                       ' - ' num2str(numel(unique(D_0))) ' intensity levels']);
                xlabel(seqSubPlots(1), 'x - direction parallel to leaf motion ');
                ylabel(seqSubPlots(1), 'z - direction perpendicular to leaf motion ');
                colorbar;
                drawnow;

                % show the leaf positions
                D_k =  sequencing.beam(i).fluence;
                for  k = 1:sequencing.beam(i).numOfShapes
                    shape_k = sequencing.beam(i).shapes(:, :, k);
                    [dimZ, dimX] = size(sequencing.beam(i).fluence);
                    seqSubPlots(4) = subplot(2, 2, 3.5, 'parent', seqFig);
                    imagesc(shape_k, 'parent', seqSubPlots(4));
                    hold(seqSubPlots(4), 'on');
                    set(seqSubPlots(4), 'YDir', 'normal');
                    xlabel(seqSubPlots(4), 'x - direction parallel to leaf motion ');
                    ylabel(seqSubPlots(4), 'z - direction perpendicular to leaf motion ');
                    title(seqSubPlots(4), ['beam # ' num2str(i) ' shape # ' num2str(k) ' d_k = ' num2str(sequencing.beam(i).shapesWeight(k))]);
                    for j = 1:dimZ
                        leftLeafIx = find(shape_k(j, :) > 0, 1, 'first');
                        rightLeafIx = find(shape_k(j, :) > 0, 1, 'last');
                        if leftLeafIx > 1
                            plot(seqSubPlots(4), [.5 leftLeafIx - .5], j - [.5 .5], 'w', 'LineWidth', 2);
                            plot(seqSubPlots(4), [.5 leftLeafIx - .5], j + [.5 .5], 'w', 'LineWidth', 2);
                            plot(seqSubPlots(4), [leftLeafIx - .5 leftLeafIx - .5], j + [.5 -.5], 'w', 'LineWidth', 2);
                        end
                        if rightLeafIx < dimX
                            plot(seqSubPlots(4), [dimX + .5 rightLeafIx + .5], j - [.5 .5], 'w', 'LineWidth', 2);
                            plot(seqSubPlots(4), [dimX + .5 rightLeafIx + .5], j + [.5 .5], 'w', 'LineWidth', 2);
                            plot(seqSubPlots(4), [rightLeafIx + .5 rightLeafIx + .5], j + [.5 -.5], 'w', 'LineWidth', 2);
                        end
                        if isempty(rightLeafIx) && isempty (leftLeafIx)
                            plot(seqSubPlots(4), [dimX + .5 .5], j - [.5 .5], 'w', 'LineWidth', 2);
                            plot(seqSubPlots(4), [dimX + .5 .5], j + [.5 .5], 'w', 'LineWidth', 2);
                            plot(seqSubPlots(4), .5 * dimX * [1 1] + 0.5, j + [.5 -.5], 'w', 'LineWidth', 2);
                        end
                    end
                    pause(1);

                    % Plot residual intensity matrix.
                    D_k = D_k - shape_k; % residual intensity matrix for visualization
                    seqSubPlots(2) = subplot(2, 2, 2, 'parent', seqFig);
                    imagesc(D_k, 'parent', seqSubPlots(2));
                    set(seqSubPlots(2), 'CLim', [0 this.numLevels], 'YDir', 'normal');
                    title(seqSubPlots(2), ['k = ' num2str(k)]);
                    colorbar;
                    drawnow;

                    axis tight;
                    drawnow;
                end

            end
        end

    end
    methods  (Static)

        function [available, msg] = isAvailable(pln, machine)

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end
            % checkBasic
            available = isfield(machine, 'meta') && isfield(machine, 'data');

            available = available && any(isfield(machine.meta, {'machine', 'name'}));

            if ~available
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
            else
                msg = [];
            end
        end

        function [pln, stf] = aperture2collimation(pln, stf, sequence, visBool)

            if nargin < 4
                visBool = false;
            end

            bixelWidth = sequence.apertureInfo.bixelWidth;
            leafWidth = bixelWidth;
            convResolution = 0.5; % [mm]

            % The collimator limits are inferred here from the apertureInfo. This could
            % be handled differently by explicitly storing collimator info in the base
            % data?
            symmetricMLClimits = vertcat(sequence.apertureInfo.beam.MLCWindow);
            symmetricMLClimits = max(abs(symmetricMLClimits));
            fieldWidth = 2 * max(symmetricMLClimits);

            % modify basic pln variables
            pln.propStf.bixelWidth = 'field';
            pln.propStf.collimation.convResolution = 0.5; % [mm]
            pln.propStf.collimation.fieldWidth = fieldWidth;
            pln.propStf.collimation.leafWidth = leafWidth;

            %
            % [bixelFieldX,bixelFieldY] = ndgrid(-fieldWidth/2:bixelWidth:fieldWidth/2,-fieldWidth/2:leafWidth:fieldWidth/2);
            [convFieldX, convFieldY] = meshgrid(-fieldWidth / 2:convResolution:fieldWidth / 2);

            % TODO: Not used in calcPhotonDose but imported from DICOM
            % pln.propStf.collimation.Devices ...
            % pln.propStf.collimation.numOfFields
            % pln.propStf.collimation.beamMeterset

            for iBeam = 1:numel(stf)
                stfTmp = stf(iBeam);
                beamSequencing = sequence.beam(iBeam);
                beamAperture = sequence.apertureInfo.beam(iBeam);

                stfTmp.bixelWidth = 'field';

                nShapes = beamSequencing.numOfShapes;

                stfTmp.numOfRays = 1; %
                stfTmp.numOfBixelsPerRay = nShapes;
                stfTmp.totalNumOfBixels = nShapes;

                ray = struct();
                ray.rayPos_bev = [0 0 0];
                ray.targetPoint_bev = [0 stfTmp.SAD 0];
                ray.weight = 1;
                ray.energy = stfTmp.ray(1).energy;
                ray.beamletCornersAtIso = stfTmp.ray(1).beamletCornersAtIso;
                ray.rayCorners_SCD = stfTmp.ray(1).rayCorners_SCD;

                % ray.shape = beamSequencing.sum;
                shapeTotalF = zeros(size(convFieldX));

                ray.shapes = struct();
                for iShape = 1:nShapes
                    currShape = beamAperture.shape(iShape);
                    activeLeafPairPosY = beamAperture.leafPairPos;
                    F = zeros(size(convFieldX));
                    if visBool
                        hF = figure;
                        imagesc(F);
                        title(sprintf('Beam %d, Shape %d', iBeam, iShape));
                        hold on;
                    end
                    for iLeafPair = 1:numel(activeLeafPairPosY)
                        posY = activeLeafPairPosY(iLeafPair);
                        ixY = convFieldY >= posY - leafWidth / 2 & convFieldY < posY + leafWidth / 2;
                        ixX = convFieldX >= currShape.leftLeafPos(iLeafPair) & convFieldX < currShape.rightLeafPos(iLeafPair);
                        ix = ixX & ixY;
                        F(ix) = 1;
                        if visBool
                            figure(hF);
                            imagesc(F);
                            drawnow;
                            pause(0.1);
                        end
                    end

                    if visBool
                        pause(1);
                        close(hF);
                    end

                    F = F * currShape.weight;
                    shapeTotalF = shapeTotalF + F;

                    ray.shapes(iShape).convFluence = F;
                    ray.shapes(iShape).shapeMap = currShape.shapeMap;
                    ray.shapes(iShape).weight = currShape.weight;
                    ray.shapes(iShape).leftLeafPos = currShape.leftLeafPos;
                    ray.shapes(iShape).rightLeafPos = currShape.rightLeafPos;
                    ray.shapes(iShape).leafPairCenterPos = activeLeafPairPosY;
                end

                ray.shape = shapeTotalF;
                ray.weight = ones(1, nShapes);
                ray.collimation = pln.propStf.collimation;
                stfTmp.ray = ray;

                stf(iBeam) = stfTmp;
            end
        end

        function [geometry, bixelIndOffset] = getMLCGeometry(stfBeam, numOfMLCLeafPairs, bixelIndOffset)
            % Derive the per-beam MLC geometry from an stf beam's ray
            % positions: bixel index map, leaf position limits, active leaf
            % pairs and MLC window. Shared by aperture-info construction in
            % the photon sequencers and the fine-angle recalculation.
            %
            % input:
            %   stfBeam:            single stf beam entry (stf(i))
            %   numOfMLCLeafPairs:  total number of physical MLC leaf pairs
            %   bixelIndOffset:     index of the last bixel of the preceding
            %                       beams; the beam's rays are numbered
            %                       starting from bixelIndOffset + 1
            %
            % output:
            %   geometry:        struct with fields numOfActiveLeafPairs,
            %                    leafPairPos, isActiveLeafPair, centralLeafPair,
            %                    limLeft, limRight, bixelIndMap, posOfCornerBixel,
            %                    MLCWindow
            %   bixelIndOffset:  updated offset (input + stfBeam.numOfRays)

            bixelWidth = stfBeam.bixelWidth; % [mm]

            % define central leaf pair (here we want the 0mm position to be in the
            % center of a leaf pair, e.g. leaf 41 stretches from -2.5mm to 2.5mm
            % for a bixel/leafWidth of 5mm and 81 leaf pairs)
            centralLeafPair = ceil(numOfMLCLeafPairs / 2);

            % get x- and z-coordinates of bixels
            rayPos_bev = reshape([stfBeam.ray.rayPos_bev], 3, []);
            X = rayPos_bev(1, :)';
            Z = rayPos_bev(3, :)';

            % create ray-map
            maxX = max(X);
            minX = min(X);
            maxZ = max(Z);
            minZ = min(Z);

            dimX = (maxX - minX) / bixelWidth + 1;
            dimZ = (maxZ - minZ) / bixelWidth + 1;

            rayMap = zeros(dimZ, dimX);

            % get indices for x and z positions
            xPos = (X - minX) / bixelWidth + 1;
            zPos = (Z - minZ) / bixelWidth + 1;

            % get indices in the ray-map
            indInRay = zPos + (xPos - 1) * dimZ;

            % fill ray-map
            rayMap(indInRay) = 1;

            % create map of bixel indices
            bixelIndMap = NaN * ones(dimZ, dimX);
            bixelIndMap(indInRay) = (1:stfBeam.numOfRays) + bixelIndOffset;
            bixelIndOffset = bixelIndOffset + stfBeam.numOfRays;

            % get leaf limits from the leaf map
            limLeft = NaN * ones(dimZ, 1);
            limRight = NaN * ones(dimZ, 1);
            for l = 1:dimZ
                lim_lInd = find(rayMap(l, :), 1, 'first');
                lim_rInd = find(rayMap(l, :), 1, 'last');
                % the physical position [mm] can be calculated from the indices
                limLeft(l) = (lim_lInd - 1) * bixelWidth + minX - 1 / 2 * bixelWidth;
                limRight(l) = (lim_rInd - 1) * bixelWidth + minX + 1 / 2 * bixelWidth;
            end

            % find upmost and downmost leaf pair
            topLeafPair = centralLeafPair - maxZ / bixelWidth;
            bottomLeafPair = centralLeafPair - minZ / bixelWidth;

            % create bool map of active leaf pairs
            isActiveLeafPair = zeros(numOfMLCLeafPairs, 1);
            isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;

            % getting the dimensions of the MLC in order to be able to plot the
            % shapes using physical coordinates
            MLCWindow = [minX - bixelWidth / 2 maxX + bixelWidth / 2 ...
                         minZ - bixelWidth / 2 maxZ + bixelWidth / 2];

            geometry.numOfActiveLeafPairs = dimZ;
            geometry.leafPairPos          = unique(Z);
            geometry.isActiveLeafPair     = isActiveLeafPair;
            geometry.centralLeafPair      = centralLeafPair;
            geometry.limLeft                = limLeft;
            geometry.limRight                = limRight;
            geometry.bixelIndMap          = bixelIndMap;
            geometry.posOfCornerBixel     = [minX 0 minZ];
            geometry.MLCWindow            = MLCWindow;
        end

        function newBeam = discardApertures(beam, numToKeep, selection)
            % The sequencing algorithm generates an a priori unknown number
            % of apertures. We only want to keep a certain number of them
            % (numToKeep), see availableApertureSelection for how they are
            % picked. The kept shapes are preserved in their original
            % trajectory order and get new weights.
            %
            % input:
            %   beam:       beam struct containing original shapes and intensities
            %   numToKeep:  number of apertures to keep
            %   selection:  'doseAreaProduct' (default) or 'leastSquares'
            %
            % output:
            %   newBeam:    beam struct with shapes and re-scaled intensities

            if nargin < 3 || isempty(selection)
                selection = 'doseAreaProduct';
            end

            numToKeep = min(numToKeep, beam.numOfShapes);

            newBeam = beam;
            newBeam.shapes = zeros(size(beam.shapes, 1), size(beam.shapes, 2), numToKeep);
            newBeam.shapesWeight = zeros(numToKeep, 1);
            newBeam.numOfShapes = numToKeep;

            if numToKeep == 0
                return
            end

            switch selection
                case 'doseAreaProduct'
                    [keepIx, weights] = matRad_PhotonLeafSequencerAbstract.selectByDoseAreaProduct(beam, numToKeep);
                case 'leastSquares'
                    [keepIx, weights] = matRad_PhotonLeafSequencerAbstract.selectByLeastSquares(beam, numToKeep);
                otherwise
                    matRad_cfg = MatRad_Config.instance();
                    matRad_cfg.dispError('Invalid aperture selection ''%s''!', selection);
            end

            % Assign into the preallocations so that the shape/class of the
            % outputs does not depend on that of the inputs (shapes stay
            % double even for a logical input, weights stay a column vector
            % even for a row-vector input).
            newBeam.shapes(:, :, :) = beam.shapes(:, :, keepIx);
            newBeam.shapesWeight(:) = weights;
        end

        function selection = availableApertureSelection()
            % Strategies for reducing the generated apertures of a beam to the
            % number that can actually be delivered::
            %
            %   doseAreaProduct  keep the apertures with the highest product of
            %                    open area and weight, and rescale them so that
            %                    the total dose-area product is maintained.
            %                    Since the stratification gives every segment
            %                    the same weight, this ranks by open area.
            %   leastSquares     keep the apertures that best reconstruct the
            %                    fluence the full decomposition would deliver,
            %                    with the weights refit to that fluence.
            selection = {'doseAreaProduct', 'leastSquares'};
        end

        function [keepIx, weights] = selectByDoseAreaProduct(beam, numToKeep)
            % Original selection: highest open-area x weight, rescaled to
            % preserve the total dose-area product.

            DAP = zeros(beam.numOfShapes, 1);
            for shape = 1:beam.numOfShapes
                DAP(shape) = nnz(beam.shapes(:, :, shape)) .* beam.shapesWeight(shape);
            end

            [~, DAPSort] = sort(DAP, 'descend');

            totDAPAll = sum(DAP(:));
            keepIx = sort(DAPSort(1:numToKeep));
            totDAPKeep = sum(DAP(keepIx));

            if totDAPKeep > 0
                weights = beam.shapesWeight(keepIx) .* (totDAPAll ./ totDAPKeep);
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning(['All %d kept apertures have zero dose-area product; ' ...
                                        'the beam is left with zero weight.'], numToKeep);
                weights = zeros(numToKeep, 1);
            end
        end

        function [keepIx, weights] = selectByLeastSquares(beam, numToKeep)
            % Pick the subset of apertures that best reproduces the fluence the
            % full decomposition delivers, and refit their weights to it.
            %
            % Forward stepwise selection: in each step try every aperture that
            % is not selected yet, refit ALL weights by non-negative least
            % squares, and keep the candidate with the smallest residual.
            % Refitting before choosing matters - picking by the residual
            % reduction of a single aperture (plain matching pursuit) commits
            % to an aperture before the weights of the others can adapt, and
            % then never recovers from overshooting.
            %
            % Cost is numToKeep x numOfShapes tiny least squares solves, which
            % is nothing next to the decomposition itself.

            numShapes = beam.numOfShapes;
            shapeMx = reshape(double(beam.shapes(:, :, 1:numShapes)), [], numShapes);

            % the fluence the beam would deliver if every aperture were kept
            % (the weights may come in as a row vector)
            allWeights = beam.shapesWeight(1:numShapes);
            targetFluence = shapeMx * allWeights(:);

            selected = zeros(numToKeep, 1);
            weights = [];

            for i = 1:numToKeep
                bestResidual = inf;
                for candidate = 1:numShapes
                    if i > 1 && any(selected(1:i - 1) == candidate)
                        continue
                    end

                    trialIx = [selected(1:i - 1); candidate];
                    trialWeights = matRad_PhotonLeafSequencerAbstract.solveNonNegativeLeastSquares( ...
                                                                                                   shapeMx(:, trialIx), targetFluence);
                    residual = norm(targetFluence - shapeMx(:, trialIx) * trialWeights);

                    if residual < bestResidual
                        bestResidual = residual;
                        selected(i) = candidate;
                        weights = trialWeights;
                    end
                end
            end

            % Least squares fits the shape of the fluence but not its integral,
            % so the kept apertures deliver systematically less than the full
            % decomposition and the plan comes out underdosed. Rescale to the
            % original dose-area product, as the heuristic does - for binary
            % aperture maps that is the same as preserving the fluence
            % integral, and it leaves the relative weights of the fit intact.
            deliveredFluence = sum(shapeMx(:, selected) * weights);
            if deliveredFluence > 0
                weights = weights * (sum(targetFluence) / deliveredFluence);
            end

            % restore the original trajectory order of the kept apertures
            [keepIx, order] = sort(selected);
            weights = weights(order);
        end

        function x = solveNonNegativeLeastSquares(designMx, rhs)
            % Minimal non-negative least squares for the handful of columns
            % used here: solve unconstrained, drop the most negative
            % coefficient, repeat. Not a full Lawson-Hanson (coefficients are
            % never re-admitted), which is adequate for numToKeep columns.

            numCols = size(designMx, 2);
            x = zeros(numCols, 1);
            active = true(numCols, 1);

            for i = 1:numCols
                ixActive = find(active);
                xActive = designMx(:, active) \ rhs;

                if all(xActive >= 0)
                    x(active) = xActive;
                    return
                end

                [~, worst] = min(xActive);
                active(ixActive(worst)) = false;

                if ~any(active)
                    return
                end
            end
        end

    end
end
