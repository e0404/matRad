function apertureInfo = matRad_leafTouching(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to improve instances of leaf touching by moving leaves
% from the centre to sweep with the non-touching leaves.
%
% Currently only works with VMAT, add option to work with IMRT (not as
% crucial)
%
% call
%   apertureInfo = matRad_leafTouching(apertureInfo)
%
% input
%   apertureInfo: matRad aperture weight and shape info struct (requires
%                 apertureInfo.arc)
%
% output
%   apertureInfo: matRad aperture weight and shape info struct
%
% Note on the two sampling modes below: the function dispatches on whether
% shape(1) already carries leftLeafPosInitial/Final, sampling the trajectory
% either once (beam centre) or twice (dose sector borders) per DAO beam. Only
% the once-per-beam mode has ever run. The sole call site
% (matRad_PhotonSequencerVMATAbstract) does not set those fields beforehand --
% this function creates them itself, at the end of its own run -- and that has
% been true since the function was first added. The twice-per-beam mode was
% presumably written for a second, post-DAO cleanup invocation that was never
% wired up; it read a field that never existed on apertureInfo.beam, so it
% would have errored on its first execution. Its output has therefore never
% been validated: do not enable a second invocation without checking the
% sweep-reset semantics (the limLeft/limRight substitution at FMO sector
% borders) against intent.
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
% Only DAO beams carry a shape at this point (interpolated beams get theirs
% from matRad_daoVec2ApertureInfo later), and beam 1 is not necessarily a DAO
% beam - with centre-sampled dose angles the arc opens with interpolated
% beams. Probe the first DAO beam instead.
refBeam = find([apertureInfo.arc.beam.isDAOBeam], 1);
dimZ = apertureInfo.beam(refBeam).numOfActiveLeafPairs;
numBeams = nnz([apertureInfo.arc.beam.isDAOBeam]);
if ~isfield(apertureInfo.beam(refBeam).shape(1), 'leftLeafPosInitial')
    % Each non-interpolated beam should have 1 left/right leaf position
    leftLeafPoss = nan(dimZ, numBeams);
    rightLeafPoss = nan(dimZ, numBeams);
    gantryAngles = zeros(1, numBeams);
else
    % Each non-interpolated beam should have 2 left/right leaf positions
    leftLeafPoss = nan(dimZ, 2 * numBeams);
    rightLeafPoss = nan(dimZ, 2 * numBeams);
    gantryAngles = zeros(1, 2 * numBeams);
end
initBorderGantryAngles = unique([apertureInfo.arc.beam.FMOAngleBorders]);
initBorderLeftLeafPoss = nan(dimZ, numel(initBorderGantryAngles));

l = 1;
m = 1;
% collect all leaf positions
for k = 1:numel(apertureInfo.beam)
    if (k ~= 1 && apertureInfo.beam(k).gantryAngle == apertureInfo.beam(k - 1).gantryAngle) || ~apertureInfo.arc.beam(k).isDAOBeam
        continue
    end

    if ~isfield(apertureInfo.beam(refBeam).shape(1), 'leftLeafPosInitial')
        leftLeafPoss(:, l) = apertureInfo.beam(k).shape(1).leftLeafPos;
        rightLeafPoss(:, l) = apertureInfo.beam(k).shape(1).rightLeafPos;
        gantryAngles(l) = apertureInfo.beam(k).gantryAngle;

        l = l + 1;
    else
        leftLeafPoss(:, l) = apertureInfo.beam(k).shape(1).leftLeafPosInitial;
        rightLeafPoss(:, l) = apertureInfo.beam(k).shape(1).rightLeafPosInitial;
        gantryAngles(l) = apertureInfo.arc.beam(k).doseAngleBorders(1);

        l = l + 1;

        leftLeafPoss(:, l) = apertureInfo.beam(k).shape(1).leftLeafPosFinal;
        rightLeafPoss(:, l) = apertureInfo.beam(k).shape(1).rightLeafPosFinal;
        gantryAngles(l) = apertureInfo.arc.beam(k).doseAngleBorders(2);

        l = l + 1;
    end

    % Only important when cleaning up instances of opposing
    % leaves touching.
    if apertureInfo.arc.beam(k).isFMOBeam
        if apertureInfo.arc.beam(k).leafDir == 1
            % This means that the current arc sector is moving
            % in the normal direction (L-R).
            initBorderLeftLeafPoss(:, m) = apertureInfo.beam(k).limLeft;

        elseif apertureInfo.arc.beam(k).leafDir == -1
            % This means that the current arc sector is moving
            % in the reverse direction (R-L).
            initBorderLeftLeafPoss(:, m) = apertureInfo.beam(k).limRight;
        end
        m = m + 1;

        % end of last sector
        if m == numel(initBorderGantryAngles)
            % This gives ending angle of the current sector.
            if apertureInfo.arc.beam(k).leafDir == 1
                % This means that the current arc sector is moving
                % in the normal direction (L-R), so the next arc
                % sector is moving opposite
                initBorderLeftLeafPoss(:, m) = apertureInfo.beam(k).limRight;
            elseif apertureInfo.arc.beam(k).leafDir == -1
                % This means that the current arc sector is moving
                % in the reverse direction (R-L), so the next
                % arc sector is moving opposite
                initBorderLeftLeafPoss(:, m) = apertureInfo.beam(k).limLeft;
            end
        end
    end
end

[gantryAngles, ind] = unique(gantryAngles);
leftLeafPoss = leftLeafPoss(:, ind);
rightLeafPoss = rightLeafPoss(:, ind);

% Any time leaf pairs are touching, they are set to
% be in the middle of the field.  Instead, move them
% so that they are still touching, but that they
% follow the motion of the MLCs across the field.
for row = 1:dimZ

    touchingInd = find(leftLeafPoss(row, :) == rightLeafPoss(row, :));

    if ~exist('leftLeafPossAug', 'var')
        % leftLeafPossAug = [reshape(mean([leftLeafPoss(:) rightLeafPoss(:)],2),size(leftLeafPoss)),borderLeftLeafPoss];
        leftLeafPossAugTemp = reshape(mean([leftLeafPoss(:) rightLeafPoss(:)], 2), size(leftLeafPoss));

        numRep = 0;
        repInd = nan(size(gantryAngles));
        for j = 1:numel(gantryAngles)
            if any(gantryAngles(j) == initBorderGantryAngles)
                % replace leaf positions with the ones at
                % the borders (eliminates repetitions)
                numRep = numRep + 1;
                % these are the gantry angles that are
                % repeated
                repInd(numRep) = j;

                delInd = find(gantryAngles(j) == initBorderGantryAngles);
                leftLeafPossAugTemp(:, j) = initBorderLeftLeafPoss(:, delInd);
                initBorderLeftLeafPoss(:, delInd) = [];
                initBorderGantryAngles(delInd) = [];
            end
        end
        repInd(isnan(repInd)) = [];
        leftLeafPossAug = [leftLeafPossAugTemp, initBorderLeftLeafPoss];
        gantryAnglesAug = [gantryAngles, initBorderGantryAngles];
    end
    % Support points are all sampled columns, not just the first numBeams:
    % the twice-per-beam mode fills two columns (initial/final) per DAO beam.
    % Equivalent to 1:numBeams in the once-per-beam mode that actually runs
    % (duplicate gantry angles are already skipped when collecting, so the
    % unique above does not shrink the array there).
    notTouchingInd = [setdiff(1:numel(gantryAngles), touchingInd), repInd];
    notTouchingInd = unique(notTouchingInd);
    % make sure to include the repeated ones in the
    % interpolation!

    notTouchingIndAug = [notTouchingInd, (1 + numel(gantryAngles)):(numel(gantryAngles) + numel(initBorderGantryAngles))];

    leftLeafPoss(row, touchingInd) = interp1(gantryAnglesAug(notTouchingIndAug), ...
                                             leftLeafPossAug(row, notTouchingIndAug), gantryAngles(touchingInd)) - 0.5;
    rightLeafPoss(row, touchingInd) = leftLeafPoss(row, touchingInd) + 1;
end

% finally, set new leaf positions
for i = 1:numel(apertureInfo.beam)
    % angles at which the leaf positions are sampled: beam centre and the
    % initial/final borders of the beam's dose angle sector
    angleC = apertureInfo.beam(i).gantryAngle;
    angleI = apertureInfo.arc.beam(i).doseAngleBorders(1);
    angleF = apertureInfo.arc.beam(i).doseAngleBorders(2);
    limL = apertureInfo.beam(i).limLeft;
    limR = apertureInfo.beam(i).limRight;

    apertureInfo.beam(i).shape(1).leftLeafPos = max((interp1(gantryAngles', leftLeafPoss', angleC))', limL);
    apertureInfo.beam(i).shape(1).rightLeafPos = min((interp1(gantryAngles', rightLeafPoss', angleC))', limR);

    apertureInfo.beam(i).shape(1).leftLeafPosInitial = max((interp1(gantryAngles', leftLeafPoss', angleI))', limL);
    apertureInfo.beam(i).shape(1).rightLeafPosInitial = min((interp1(gantryAngles', rightLeafPoss', angleI))', limR);

    apertureInfo.beam(i).shape(1).leftLeafPosFinal = max((interp1(gantryAngles', leftLeafPoss', angleF))', limL);
    apertureInfo.beam(i).shape(1).rightLeafPosFinal = min((interp1(gantryAngles', rightLeafPoss', angleF))', limR);
end

end
