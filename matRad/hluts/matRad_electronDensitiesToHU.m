function ct = matRad_electronDensitiesToHU(ct, radiationMode)
% matRad function to calculate recalculate HU values from equivalent densities
% This is done to provide downward compatibility to previous matRad
% versions where HU values were not automatically saved during the
% import process. HU values can only be calculated if the densities of the
% HLUT do not decrease with increasing HU. If several HU values map to the
% same density (e.g. the air plateau of the default HLUT), the lowest of
% these HU values is used and a warning is issued.
%
% call:
%   ct = matRad_electronDensitiesToHU(ct)
%   ct = matRad_electronDensitiesToHU(ct, radiationMode)
%
% input:
%   ct:             matRad ct struct containing cube and all additional
%                   information. If it contains an hlut (ct.hlut), this is
%                   used and radiationMode is ignored.
%   radiationMode:  optional, radiationMode as character array
%                   (e.g. 'photons') or pln-struct with the field
%                   radiationMode, used to select the HLUT if ct.hlut does
%                   not exist. Defaults to 'photons', where the HLUT maps
%                   HU to relative electron density. For particles, the
%                   HLUT maps HU to relative stopping power instead, so the
%                   radiation mode needs to be passed explicitly.
%
% output:
%   ct: ct struct with HU and equivalent density cube
%
% References
%   -
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

matRad_cfg = MatRad_Config.instance();

if nargin < 2 || isempty(radiationMode)
    radiationMode = 'photons';
end

% load hlut
if isfield(ct, 'hlut') % if hlut stored upon import use this one!
    hlut = ct.hlut;
else
    hlut = matRad_loadHLUT(ct, radiationMode);
end

% the HLUT can only be inverted if the densities do not decrease with HU
hlut = sortrows(hlut, [1 2]);
if any(diff(hlut(:, 2)) < 0)
    matRad_cfg.dispWarning('Reconversion of HU values could not be done because the densities of the HLUT decrease with HU.');
    return
end

% plateaus (several HU values with the same density) are not bijective.
% Densities on a plateau are mapped to its lowest HU value, densities
% between two HLUT points are interpolated from the highest HU value of the
% lower point to the lowest HU value of the upper point, so that the
% segments adjacent to a plateau are inverted exactly. 'first' / 'last' are
% needed for Octave.
[densities, firstIx] = unique(hlut(:, 2), 'first');
[~, lastIx] = unique(hlut(:, 2), 'last');
huLow = hlut(firstIx, 1);
huHigh = hlut(lastIx, 1);
if numel(densities) < size(hlut, 1)
    matRad_cfg.dispWarning(['HLUT is not bijective, several HU values map to the same density. ', ...
                            'Using the lowest HU value for each of these densities.']);
end

for i = 1:ct.numOfCtScen
    % Manual adjustments if ct data is corrupt. If some values are out of range
    % of the LUT, then these values are adjusted.
    if max(ct.cube{i}(:)) > densities(end)
        matRad_cfg.dispWarning('projecting out of range electron density values');
        ct.cube{i}(ct.cube{i}(:) > densities(end)) = densities(end);
    end
    if min(ct.cube{i}(:)) < densities(1)
        matRad_cfg.dispWarning('projecting out of range electron density values');
        ct.cube{i}(ct.cube{i}(:) < densities(1)) = densities(1);
    end

    ct.cubeHU{i} = matRad_invertHlut(densities, huLow, huHigh, ct.cube{i});
end

end

function hu = matRad_invertHlut(densities, huLow, huHigh, cube)
% piecewise linear inverse of an HLUT with plateaus, see above. cube is
% expected to lie within [densities(1), densities(end)]
hu = NaN(size(cube));

% densities on an HLUT point (including plateaus) get its lowest HU value
[isPoint, pointIx] = ismember(cube, densities);
hu(isPoint) = huLow(pointIx(isPoint));

% all other densities lie strictly within a segment [k, k+1], NaN stays NaN
isInSegment = ~isPoint & ~isnan(cube);
if any(isInSegment(:))
    q = cube(isInSegment);
    k = floor(interp1(densities, (1:numel(densities))', q));
    k = min(k, numel(densities) - 1);
    t = (q - densities(k)) ./ (densities(k + 1) - densities(k));
    hu(isInSegment) = huHigh(k) + t .* (huLow(k + 1) - huHigh(k));
end

end
