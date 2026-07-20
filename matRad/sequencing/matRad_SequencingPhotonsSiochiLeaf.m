classdef  matRad_SequencingPhotonsSiochiLeaf < matRad_PhotonSequencerVMATAbstract

    % matRad_SequencingPhotonsSiochiLeaf: photon MLC leaf sequencing after
    %   Siochi (1999), using rod pushing with collision and tongue-and-groove
    %   correction to convert the optimized fluence into deliverable segments.
    %   Supports both static (IMRT) and VMAT (dynamic/arc) delivery - see
    %   this.runVMAT and matRad_PhotonSequencerVMATAbstract.
    %
    % References
    %   [1] https://www.ncbi.nlm.nih.gov/pubmed/10078655

    properties (Constant)
        name = 'Photons Siochi Leaf Sequencer'
        shortName = 'siochi'
        possibleRadiationModes = {'photons'}

    end

    methods

        function this = matRad_SequencingPhotonsSiochiLeaf(pln)
            if nargin < 1
                pln = [];
            end
            this = this@matRad_PhotonSequencerVMATAbstract(pln);
        end

        function sequence = sequence(this, w, stf)
            if this.runVMAT
                sequence = this.sequenceDynamic(w, stf);
            else
                sequence = this.sequenceStatic(w, stf);
            end

            if this.preconditioner
                sequence.apertureInfo = matRad_preconditionFactors(sequence.apertureInfo);
            end

            if this.visMode
                this.plotSegments(sequence);
            end
        end

        function sequence = sequenceStatic(this, w, stf)

            offset = 0;

            for i = 1:numel(stf)

                [D_0, D_k, shapes, calFac, indInMx] = this.initBeam(stf(i), w(1 + offset:stf(i).numOfRays + offset));

                shapesWeight = zeros(10000, 1);
                k = 0;
                if sum(w(1 + offset:stf(i).numOfRays + offset)) > 0

                    % Decompose the port, do rod pushing
                    [tops, bases] = this.decomposePort(D_k);
                    % Form segments
                    [shapes, shapesWeight, k] = this.convertToSegments(shapes, shapesWeight, k, tops, bases);

                    sequence.beam(i).numOfShapes  = k;
                    sequence.beam(i).shapes       = shapes(:, :, 1:k);
                    sequence.beam(i).shapesWeight = shapesWeight(1:k) / this.sequencingLevel * calFac;
                    sequence.beam(i).bixelIx      = 1 + offset:size(stf(i).ray, 2) + offset;
                    sequence.beam(i).fluence      = D_0;
                    sequence.beam(i).sum          = zeros(size(D_0));

                    for j = 1:k
                        sequence.beam(i).sum = sequence.beam(i).sum + sequence.beam(i).shapes(:, :, j) * sequence.beam(i).shapesWeight(j);
                    end
                else
                    sequence.beam(i).numOfShapes  = 1;
                    sequence.beam(i).shapes       = zeros(size(D_0));
                    sequence.beam(i).shapesWeight = zeros(size(D_0));
                    sequence.beam(i).bixelIx      = 1 + offset:size(stf(i).ray, 2) + offset;
                    sequence.beam(i).fluence      = zeros(size(D_0));
                    sequence.beam(i).sum          = zeros(size(D_0));
                end
                if stf(i).numOfRays > 1
                    sequence.w(1 + offset:stf(i).numOfRays + offset, 1) = sequence.beam(i).sum(indInMx);
                else
                    sequence.w(1 + offset:stf(i).numOfRays + offset, 1) = w(1 + offset:stf(i).numOfRays + offset);
                end
                offset = offset + stf(i).numOfRays;

            end

            sequence = this.sequencing2ApertureInfo(sequence, stf);
        end

        function sequence = sequenceDynamic(this, w, stf)
            % VMAT (dynamic/arc) sequencing: gates to FMO-anchor beams,
            % smooths the fluence before decomposition, keeps re-decomposing
            % at increasing sequencingLevel until enough shapes exist for
            % this beam's DAO-angle children, then spreads the result across
            % those children and builds the VMAT apertureInfo. Ported from
            % the former matRad_siochiLeafSequencing.m functional
            % implementation.

            offset = 0;

            for i = 1:numel(stf)
                numOfRaysPerBeam = stf(i).numOfRays;

                if ~stf(i).arc.isFMOBeam
                    sequence.w(1 + offset:numOfRaysPerBeam + offset, 1) = 0;
                    sequence.beam(i).bixelIx = 1 + offset:numOfRaysPerBeam + offset;
                    offset = offset + numOfRaysPerBeam;
                    continue % this beam carries no shapes of its own; arcSequencing fills it in as a child of an FMO beam
                end
                numToKeep = stf(i).arc.numOfChildren;

                wOfCurrBeams = w(1 + offset:numOfRaysPerBeam + offset) .* ones(size(stf(i).ray, 2), 1);

                X = ones(size(stf(i).ray, 2), 1) * NaN;
                Z = ones(size(stf(i).ray, 2), 1) * NaN;
                for j = 1:size(stf(i).ray, 2)
                    X(j) = stf(i).ray(j).rayPos_bev(:, 1);
                    Z(j) = stf(i).ray(j).rayPos_bev(:, 3);
                end

                minX = min(X);
                maxX = max(X);
                minZ = min(Z);
                maxZ = max(Z);

                dimOfFluenceMxX = (maxX - minX) / stf(i).bixelWidth + 1;
                dimOfFluenceMxZ = (maxZ - minZ) / stf(i).bixelWidth + 1;

                fluenceMx = zeros(dimOfFluenceMxZ, dimOfFluenceMxX);

                xPos = (X - minX) / stf(i).bixelWidth + 1;
                zPos = (Z - minZ) / stf(i).bixelWidth + 1;

                indInFluenceMx = zPos + (xPos - 1) * dimOfFluenceMxZ;
                fluenceMx(indInFluenceMx) = wOfCurrBeams;

                % Gaussian fluence smoothing (VMAT only - would corrupt static
                % IMRT fluence reproduction otherwise)
                fluenceMx = this.smoothFluenceForArc(fluenceMx);

                numOfLevels = this.sequencingLevel;
                notFinished = true;

                if sum(wOfCurrBeams) > 0
                    while notFinished
                        % keep re-decomposing at an increasing number of
                        % levels until there are at least as many shapes as
                        % this FMO beam's DAO-angle children (numToKeep)

                        calFac = max(fluenceMx(:));
                        D_k = round(fluenceMx / calFac * numOfLevels);
                        D_0 = D_k;

                        shapes = NaN * ones(dimOfFluenceMxZ, dimOfFluenceMxX, 10000);
                        shapesWeight = zeros(10000, 1);
                        k = 0;

                        [tops, bases] = this.decomposePort(D_k);
                        [shapes, shapesWeight, k] = this.convertToSegments(shapes, shapesWeight, k, tops, bases);

                        if numToKeep ~= 0 && k < numToKeep
                            numOfLevels = numOfLevels + 1;
                        else
                            notFinished = false;
                        end
                    end

                    sequence.beam(i).numOfShapes  = k;
                    sequence.beam(i).shapes       = shapes(:, :, 1:k);
                    sequence.beam(i).shapesWeight = shapesWeight(1:k) / numOfLevels * calFac;
                    sequence.beam(i).bixelIx      = 1 + offset:numOfRaysPerBeam + offset;
                    sequence.beam(i).fluence      = D_0;
                else
                    sequence.beam(i).numOfShapes  = 1;
                    sequence.beam(i).shapes       = zeros(dimOfFluenceMxZ, dimOfFluenceMxX);
                    sequence.beam(i).shapesWeight = zeros(dimOfFluenceMxZ, dimOfFluenceMxX);
                    sequence.beam(i).bixelIx      = 1 + offset:numOfRaysPerBeam + offset;
                    sequence.beam(i).fluence      = zeros(dimOfFluenceMxZ, dimOfFluenceMxX);
                end

                if numToKeep ~= 0
                    sequence.beam(i) = this.discardApertures(sequence.beam(i), numToKeep);
                end

                sequence.beam(i).sum = zeros(dimOfFluenceMxZ, dimOfFluenceMxX);
                for shape = 1:sequence.beam(i).numOfShapes
                    sequence.beam(i).sum = sequence.beam(i).sum + sequence.beam(i).shapes(:, :, shape) * sequence.beam(i).shapesWeight(shape);
                end

                offset = offset + numOfRaysPerBeam;
            end

            % spread shapes to DAO-angle children, compute gantry rotation/MU rate
            sequence = this.applyArcSequencing(sequence, stf);

            % build the full VMAT apertureInfo and run the post-processing pipeline
            sequence.apertureInfo = this.buildVMATApertureInfo(sequence, stf);
            sequence.apertureInfo = this.postProcessVMATApertureInfo(sequence.apertureInfo);

            sequence.w = sequence.apertureInfo.bixelWeights;
        end

        function [tops, bases] = decomposePort(this, map)
            % Returns tops and bases of a fluence matrix "map" for Siochi leaf
            % sequencing algorithm (rod pushing part).  Accounts for collisions and
            % tongue and groove (Tng) effects.

            [dimZ, dimX] = size(map);
            map_nonZero = (map ~= 0);

            [D_k_Z, D_k_X] = ind2sub([dimZ, dimX], find(map_nonZero));
            minZ = min(D_k_Z);
            maxZ = max(D_k_Z);
            minX = min(D_k_X);
            maxX = max(D_k_X);

            tops = zeros(dimZ, dimX);
            bases = zeros(dimZ, dimX);

            for i = minX:maxX
                maxTop = -1;
                TnG = 1;
                for j = minZ:maxZ
                    if i == minX
                        bases(j, i) = 1;
                        tops(j, i) = bases(j, i) + map(j, i) - 1;
                    else % assign trial base positions
                        if map(j, i) >= map(j, i - 1) % current rod >= previous, match the bases
                            bases(j, i) = bases(j, i - 1);
                            tops(j, i) = bases(j, i) + map(j, i) - 1;
                        else % current rod <previous
                            if map(j, i) == 0 % rod length=0, put in in next slab after top of previous
                                bases(j, i) = tops(j, i - 1) + 1;
                                tops(j, i) = bases(j, i) - 1;
                            else % rod length~=0, match tops
                                tops(j, i) = tops(j, i - 1);
                                bases(j, i) = tops(j, i) - map(j, i) + 1;
                            end
                        end
                    end
                    % determine which rod has the highest top in column
                    if tops(j, i) > maxTop
                        maxTop = tops(j, i);
                        maxRow = j;
                    end
                end

                % Correct for collision and tongue and groove error
                while TnG
                    % go from maxRow down checking for TnG.  This occurs when a shorter
                    % rod is "peeking over" a longer one in the direction transverse to
                    % the leaf motion.  To fix this, match either the tops or bases of
                    % the rods.
                    [tops, bases] = this.fixTnGPass(tops, bases, map, i, (maxRow - 1):-1:minZ, 1);
                    % go from maxRow up checking for TnG
                    [tops, bases] = this.fixTnGPass(tops, bases, map, i, (maxRow + 1):maxZ, -1);
                    % now check if all TnG conditions have been removed
                    TnG = 0;
                    for j = (minZ + 1):maxZ
                        if map(j, i) < map(j - 1, i)
                            if tops(j, i) > tops(j - 1, i)
                                TnG = 1;
                            elseif bases(j, i) < bases(j - 1, i)
                                TnG = 1;
                            end
                        else
                            if tops(j, i) < tops(j - 1, i)
                                TnG = 1;
                            elseif bases(j, i) > bases(j - 1, i)
                                TnG = 1;
                            end
                        end
                    end
                end
            end
        end

        function [tops, bases] = fixTnGPass(~, tops, bases, map, i, jRange, neighborOffset)
            % One directional tongue-and-groove correction pass over jRange,
            % matching each rod j against its neighbor j+neighborOffset.
            % Shared by decomposePort's "down" (neighborOffset=1, comparing
            % against j+1) and "up" (neighborOffset=-1, comparing against
            % j-1) passes, which are otherwise identical.

            for j = jRange
                jn = j + neighborOffset;
                if map(j, i) < map(jn, i)
                    if tops(j, i) > tops(jn, i)
                        tops(jn, i) = tops(j, i);
                        bases(jn, i) = tops(jn, i) - map(jn, i) + 1;
                    elseif bases(j, i) < bases(jn, i)
                        bases(j, i) = bases(jn, i);
                        tops(j, i) = bases(j, i) + map(j, i) - 1;
                    end
                else
                    if tops(j, i) < tops(jn, i)
                        tops(j, i) = tops(jn, i);
                        bases(j, i) = tops(j, i) - map(j, i) + 1;
                    elseif bases(j, i) > bases(jn, i)
                        bases(jn, i) = bases(j, i);
                        tops(jn, i) = bases(jn, i) + map(jn, i) - 1;
                    end
                end
            end
        end

        function [shapes, shapesWeight, k] = convertToSegments(this, shapes, shapesWeight, k, tops, bases)
            % Convert tops and bases to shape matrices.  These are taken as to be the
            % shapes of uniform level/elevation after the rods are pushed.

            levels = max(tops(:));

            for level = 1:levels
                % check if slab is new
                if this.differentSlab(tops, bases, level)
                    k = k + 1; % increment number of unique slabs
                    shape_k = (bases <= level) .* (level <= tops); % shape of current slab
                    shapes(:, :, k) = shape_k;
                end
                shapesWeight(k) = shapesWeight(k) + 1; % if slab is not unique, this increments weight again
            end
        end

        function diffSlab = differentSlab(~, tops, bases, level)

            % Returns 1 if slab level is different than slab level-1 0 otherwise

            if level == 1 % first slab is automatically different
                diffSlab = 1;
            else
                shapeLevel = (bases <= level) .* (level <= tops); % shape of slab with current level
                shapeLevel_1 = (bases <= level - 1) .* (level - 1 <= tops); % shape of slab with previous level
                diffSlab = ~isequal(shapeLevel, shapeLevel_1); % tests if slabs are equal; isequaln was not giving correct results
            end
        end

    end
    methods  (Static)

        function [available, msg] = isAvailable(pln, machine)
            % see superclass for information

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available, msg] = matRad_PhotonSequencerAbstract.isAvailable(pln, machine);

            if ~available
                return
            else
                available = false;
                msg = [];
            end

            % checkBasic
            try
                checkBasic = isfield(machine, 'meta') && isfield(machine, 'data');

                % check modality
                checkModality = any(strcmp(matRad_SequencingPhotonsSiochiLeaf.possibleRadiationModes, machine.meta.radiationMode)) && ...
                    any(strcmp(matRad_SequencingPhotonsSiochiLeaf.possibleRadiationModes, pln.radiationMode));

                % Sanity check compatibility
                if checkModality
                    checkModality = strcmp(machine.meta.radiationMode, pln.radiationMode);
                end

                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return
            end

            available = preCheck;
        end

    end
end
