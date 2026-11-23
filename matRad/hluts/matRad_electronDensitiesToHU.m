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
isInvertible = all(diff(hlut(:, 2)) >= 0);

% plateaus (several HU values with the same density) are not bijective:
% use the lowest HU value of each plateau ('first' is needed for Octave)
if isInvertible
    [~, uniqueIx] = unique(hlut(:, 2), 'first');
    if numel(uniqueIx) < size(hlut, 1)
        matRad_cfg.dispWarning(['HLUT is not bijective, several HU values map to the same density. ', ...
                                'Using the lowest HU value for each of these densities.']);
        hlut = hlut(uniqueIx, :);
    end
end

% interpolate rel. electron dens. to HU based on lookup table
if isInvertible

    for i = 1:ct.numOfCtScen
        % Manual adjustments if ct data is corrupt. If some values are out of range
        % of the LUT, then these values are adjusted.
        if max(ct.cube{i}(:)) > max(hlut(:, 2))
            matRad_cfg.dispWarning('projecting out of range electron density values');
            ct.cube{i}(ct.cube{i}(:) > max(hlut(:, 2))) = max(hlut(:, 2));
        end
        if min(ct.cube{i}(:)) < min(hlut(:, 2))
            matRad_cfg.dispWarning('projecting out of range electron density values');
            ct.cube{i}(ct.cube{i}(:) < min(hlut(:, 2))) = min(hlut(:, 2));
        end

        ct.cubeHU{i} = interp1(hlut(:, 2), hlut(:, 1), ct.cube{i});
    end

else
    matRad_cfg.dispWarning('Reconversion of HU values could not be done because the densities of the HLUT decrease with HU.');
end
