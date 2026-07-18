function [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
% matRad function to generate a vector representation of the aperture
% weights and shapes and (optional) some meta information needed during
% optimization
%
% call
%   [apertureInfoVec, mappingMx, limMx] = matRad_daoApertureInfo2Vec(apertureInfo)
%
% input
%   apertureInfo:    aperture weight and shape info struct
%
% output
%   apertureInfoVec: vector representation of the apertue weights and shapes
%   mappingMx:       mapping of vector components to beams, shapes and leaves
%   limMx:           bounds on vector components, i.e., minimum and maximum
%                    aperture weights (0/inf) and leav positions (custom)
%
% References
%
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

% function to create a single vector for the direct aperture optimization
% first: aperture weights
% second: left leaf positions
% third: right leaf positions
% fourth: times between successive DAO gantry angles

% This is the VMAT version, with gantry time entries appended to the
% vector. Static (step-and-shoot) plans must use
% matRad_OptimizationProblemDAO's version.
if ~apertureInfo.runVMAT
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError(['Static apertureInfo passed to the VMAT vector conversion! ' ...
                          'Use matRad_OptimizationProblemDAO.matRad_daoApertureInfo2Vec (or instance dispatch) instead.']);
end

% initializing variables

% extra set of (apertureInfo.totalNumOfShapes) elements, allowing arc
% sector times to be optimized
vecLength = (apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2) + apertureInfo.totalNumOfShapes;

apertureInfoVec = NaN * ones(vecLength, 1);

offset = 0;

%% 1. aperture weights
for i = 1:size(apertureInfo.beam, 2)
    for j = 1:apertureInfo.beam(i).numOfShapes

        % In VMAT, this weight is "spread" over unoptimized beams (assume
        % constant dose rate over sector)
        apertureInfoVec(offset + j) = apertureInfo.beam(i).shape(j).jacobiScale * apertureInfo.beam(i).shape(j).weight;

    end
    offset = offset + apertureInfo.beam(i).numOfShapes;
end

% 2. left and right leaf positions
%% fill the vector for all shapes of all beams
for i = 1:size(apertureInfo.beam, 2)
    for j = 1:apertureInfo.beam(i).numOfShapes

        leafPairIx = 1:apertureInfo.beam(i).numOfActiveLeafPairs;

        if ~apertureInfo.continuousAperture

            apertureInfoVec(offset + leafPairIx) = apertureInfo.beam(i).shape(j).leftLeafPos;
            apertureInfoVec(offset + leafPairIx + apertureInfo.totalNumOfLeafPairs) = apertureInfo.beam(i).shape(j).rightLeafPos;

            offset = offset + apertureInfo.beam(i).numOfActiveLeafPairs;
        else

            if apertureInfo.propVMAT.beam(i).doseAngleDAO(1)
                apertureInfoVec(offset + leafPairIx) = apertureInfo.beam(i).shape(j).leftLeafPos_I;
                apertureInfoVec(offset + leafPairIx + apertureInfo.totalNumOfLeafPairs) = apertureInfo.beam(i).shape(j).rightLeafPos_I;

                offset = offset + apertureInfo.beam(i).numOfActiveLeafPairs;
            end

            if apertureInfo.propVMAT.beam(i).doseAngleDAO(2)
                apertureInfoVec(offset + leafPairIx) = apertureInfo.beam(i).shape(j).leftLeafPos_F;
                apertureInfoVec(offset + leafPairIx + apertureInfo.totalNumOfLeafPairs) = apertureInfo.beam(i).shape(j).rightLeafPos_F;

                offset = offset + apertureInfo.beam(i).numOfActiveLeafPairs;
            end
        end

    end
end
%% 3. time of arc sector/beam
offset = offset + apertureInfo.totalNumOfLeafPairs;

% this gives a vector of the arc lengths belonging to each optimized CP
% unique gets rid of double-counted angles (which is every interior
% angle)

optInd = [apertureInfo.propVMAT.beam.DAOBeam];
optAngleLengths = [apertureInfo.propVMAT.beam(optInd).DAOAngleBordersDiff];
optGantryRot = [apertureInfo.beam(optInd).gantryRot];
apertureInfoVec((offset + 1):end) = optAngleLengths ./ optGantryRot; % entries are the times until the next opt gantry angle is reached

%% 4. create additional information for later use
if nargout > 1

    mappingMx = NaN * ones(vecLength, 4);
    limMx = NaN * ones(vecLength, 2);

    limMx(1:(apertureInfo.totalNumOfShapes), :) = ones(apertureInfo.totalNumOfShapes, 1) * [0 inf];

    counter = 1;

    for i = 1:numel(apertureInfo.beam)
        for j = 1:apertureInfo.beam(i).numOfShapes
            mappingMx(counter, 1) = i;

            % minimum/maximum time interval between two optimized beams/gantry angles
            timeLimL = diff(apertureInfo.propVMAT.beam(i).DAOAngleBorders) / apertureInfo.propVMAT.constraints.gantryRotationSpeed(2);
            timeLimU = diff(apertureInfo.propVMAT.beam(i).DAOAngleBorders) / apertureInfo.propVMAT.constraints.gantryRotationSpeed(1);

            mappingMx(counter + (apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2), 1) = i;
            limMx(counter + (apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2), :) = [timeLimL timeLimU];

            counter = counter + 1;
        end
    end

    shapeOffset = 0;
    for i = 1:numel(apertureInfo.beam)
        for j = 1:apertureInfo.beam(i).numOfShapes
            for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs
                mappingMx(counter, 1) = i;
                mappingMx(counter, 2) = j + shapeOffset; % store global shape number for grad calc
                mappingMx(counter, 3) = j; % store local shape number
                mappingMx(counter, 4) = k; % store local leaf number

                limMx(counter, 1)     = apertureInfo.beam(i).lim_l(k);
                limMx(counter, 2)     = apertureInfo.beam(i).lim_r(k);
                counter = counter + 1;

                if apertureInfo.continuousAperture && nnz(apertureInfo.propVMAT.beam(i).doseAngleDAO) == 2
                    % redo for initial and final leaf positions
                    % might have to revisit this after looking at gradient,
                    % esp. mappingMx(counter,2)
                    % only an issue for non-interpolated deliveries
                    mappingMx(counter, 1) = i;
                    mappingMx(counter, 2) = j + shapeOffset; % store global shape number for grad calc
                    mappingMx(counter, 3) = j; % store local shape number
                    mappingMx(counter, 4) = k; % store local leaf number

                    limMx(counter, 1)     = apertureInfo.beam(i).lim_l(k);
                    limMx(counter, 2)     = apertureInfo.beam(i).lim_r(k);
                    counter = counter + 1;
                end
            end
        end
        shapeOffset = shapeOffset + apertureInfo.beam(i).numOfShapes;
    end

    lastRow = apertureInfo.totalNumOfShapes + apertureInfo.totalNumOfLeafPairs * 2;
    mappingMx(counter:lastRow, :) = mappingMx(apertureInfo.totalNumOfShapes + 1:counter - 1, :);
    limMx(counter:lastRow, :)     = limMx(apertureInfo.totalNumOfShapes + 1:counter - 1, :);

end

end
