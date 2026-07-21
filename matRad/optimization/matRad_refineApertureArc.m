function [stf, apertureInfo] = matRad_refineApertureArc(ct, cst, pln, apertureInfo, angularSpacing, varargin)
% matRad_refineApertureArc resamples an optimized VMAT plan onto a new
%   (usually finer) gantry-angle grid. The apertures optimized at the DAO
%   control points are interpolated onto the new angles and the bixel
%   weights are recomputed, yielding a steering struct and aperture info
%   ready for a forward dose calculation. No dose is computed here.
%
% call:
%   [stf, apertureInfo] = matRad_refineApertureArc(ct, cst, pln, apertureInfo, angularSpacing)
%   [stf, apertureInfo] = matRad_refineApertureArc(..., 'interpNew', tf)
%   [stf, apertureInfo] = matRad_refineApertureArc(..., 'reuseDij', tf)
%   [stf, apertureInfo] = matRad_refineApertureArc(..., 'continuousAperture', tf)
%
% input:
%   ct:              matRad ct struct
%   cst:             matRad cst cell array
%   pln:             the plan used for the original optimization (its
%                    propStf.gantryAngles hold the arc anchor points)
%   apertureInfo:    aperture info from the optimization result
%   angularSpacing:  target gantry angle spacing of the new grid [deg]
%
% name-value parameters:
%   interpNew:           interpolate new apertures at the new angles (true,
%                        default) or force every new angle to be its own DAO
%                        control point (false)
%   reuseDij:            duplicate/redirect the new beams onto the original
%                        DAO angles so that an existing dij's columns can be
%                        recycled instead of recomputing the dose. Default
%                        false. See the VMAT example for the back-projection.
%   continuousAperture:  interpolate leaf positions within a Dij arc. Default
%                        taken from apertureInfo.continuousAperture.
%
% output:
%   stf:             steering struct on the new angular grid
%   apertureInfo:    aperture info resampled onto the new grid, including the
%                    recomputed bixelWeights
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser();
p.addRequired('ct', @isstruct);
p.addRequired('cst', @iscell);
p.addRequired('pln', @isstruct);
p.addRequired('apertureInfo', @isstruct);
p.addRequired('angularSpacing', @isnumeric);
p.addParameter('interpNew', true, @islogical);
p.addParameter('reuseDij', false, @islogical);
p.addParameter('continuousAperture', ...
               matRad_getFieldOrDefault(apertureInfo, 'continuousAperture', false), @islogical);
p.parse(ct, cst, pln, apertureInfo, angularSpacing, varargin{:});

interpNew = p.Results.interpNew;
reuseDij = p.Results.reuseDij;
continuousAperture = p.Results.continuousAperture;

% Old fine angles are the DAO beams from the original optimization. They
% live in apertureInfo, so we do not need pln.propStf for this.
oldFineAngles = [apertureInfo.beam.gantryAngle];

% Set the resampling resolution on a copy of the optimization plan. The arc
% anchor points stay in pln.propStf.gantryAngles; the fine angle grid is
% computed internally by matRad_StfGeneratorPhotonVMAT.
pln.propStf.maxGantryAngleSpacing = angularSpacing;

% If not interpolating new apertures, force every new angle to be a DAO
% control point by equalising the two spacing parameters.
if ~interpNew
    pln.propStf.maxDAOGantryAngleSpacing = angularSpacing;
end

stf = matRad_generateStf(ct, cst, pln);

% -----------------------------------------------------------------------
% Dij reuse: when an existing dij is to be recycled (reuseDij) or every new
% angle is a DAO point (~interpNew), redirect the new beams onto the old DAO
% angles so that the existing Dij columns line up with the new apertures.
% -----------------------------------------------------------------------
if ~interpNew || reuseDij

    newFineAngles = [stf.gantryAngle];

    % An angle that is exactly equidistant between two old DAO angles needs
    % two stf entries, one per side, so both neighbouring Dij columns can be
    % used.
    duplicate = false(size(newFineAngles));
    for i = 1:numel(newFineAngles)
        diffs = abs(newFineAngles(i) - oldFineAngles);
        duplicate(i) = sum(diffs == min(diffs)) > 1;
    end

    tempStf = stf;
    [tempStf(:).copyInd] = deal([]);
    [tempStf(:).stfCorr] = deal([]);

    j = 1;
    for i = 1:numel(newFineAngles)
        if duplicate(i)
            % Left copy: pretend this beam sits at the previous beam angle
            tempStf(j) = stf(i);
            tempStf(j).gantryAngle = stf(i - 1).gantryAngle;
            tempStf(j).copyInd = 1;
            tempStf(j).stfCorr = false;
            j = j + 1;

            % Right copy: pretend this beam sits at the next beam angle
            tempStf(j) = stf(i);
            tempStf(j).gantryAngle = stf(i + 1).gantryAngle;
            tempStf(j).copyInd = 2;
            tempStf(j).stfCorr = false;
        else
            tempStf(j) = stf(i);
            tempStf(j).copyInd = [];
            tempStf(j).stfCorr = true;
        end
        j = j + 1;
    end
    stf = tempStf;

    % Redirect each new beam to the nearest old beam so an existing Dij
    % column can be recycled without recomputation.
    tempStf = stf;
    for i = 1:numel(stf)
        diffs = abs(stf(i).gantryAngle - oldFineAngles);
        minDiff = min(diffs);
        nearIdx = find(diffs == minDiff);   % 1 or 2 indices into oldFineAngles

        % Find where those old angles appear in the (post-duplicate) stf
        newAngles = [tempStf.gantryAngle];
        minInd1 = find(newAngles == oldFineAngles(nearIdx(1)), 1);
        minInd2 = find(newAngles == oldFineAngles(nearIdx(end)), 1);

        if reuseDij
            % Replace this beam with the nearest existing beam entirely
            if isempty(stf(i).copyInd) || stf(i).copyInd == 1
                stf(i) = tempStf(minInd1);
            elseif stf(i).copyInd == 2
                stf(i) = tempStf(minInd2);
            end
        elseif ~interpNew
            % Keep the beam but correct its angle if it is equidistant
            if numel(nearIdx) > 1
                stf(i).gantryAngle = tempStf(i).gantryAngle;
            end
        end
    end
end

% Rebuild the aperture info on the new grid and recompute the bixel weights.
% matRad_recalcApertureInfo consumes a small options struct.
recalc.pln = pln;
recalc.stf = stf;
recalc.interpNew = interpNew;
recalc = matRad_recalcApertureInfo(recalc, apertureInfo);

recalc.apertureInfo.continuousAperture = continuousAperture;
apertureInfo = matRad_recalcApertureBixelWeights(recalc.apertureInfo);

stf = recalc.stf;

end
