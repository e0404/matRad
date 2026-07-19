classdef  (Abstract) matRad_PhotonSequencerAbstract < matRad_SequencerBase

    % matRad_PhotonSequencerAbstract: Abstract base class for photon multileaf
    %   collimator (MLC) leaf sequencers. Provides the shared fluence
    %   stratification, aperture-info generation and segment visualization;
    %   concrete subclasses implement the specific leaf sequencing algorithm.
    properties
        numOfMLCLeafPairs = 80
        sequencingLevel = 5
        preconditioner = false % apply matRad_preconditionFactors to the resulting apertureInfo
    end

    methods

        function this = matRad_PhotonSequencerAbstract(pln)
            % Constructor, forwards to the base class. An explicit constructor
            % chain is required so that property assignment from pln persists
            % under Octave.
            if nargin < 1
                pln = [];
            end
            this = this@matRad_SequencerBase(pln);
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
            dCurrent = round(fluenceMx / calFac * this.sequencingLevel);

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
                [geometry, bixelIndOffset] = matRad_getMLCGeometry(stf(i), this.numOfMLCLeafPairs, bixelIndOffset);
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
                            leftLeafPos(l) = (geometry.lim_l(l) + geometry.lim_r(l)) / 2;
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
                sequence.apertureInfo.beam(i).lim_l = geometry.lim_l;
                sequence.apertureInfo.beam(i).lim_r = geometry.lim_r;
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
                set(seqSubPlots(1), 'CLim', [0 this.sequencingLevel], 'YDir', 'normal');
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
                    set(seqSubPlots(2), 'CLim', [0 this.sequencingLevel], 'YDir', 'normal');
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

    end
end
