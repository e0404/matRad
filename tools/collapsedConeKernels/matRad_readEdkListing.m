function edk = matRad_readEdkListing(fileName)
% matRad_readEdkListing reads the conical energy deposition kernel table
% from an EGSnrc EDKnrc listing file (.egslst)
%
% call:
%   edk = matRad_readEdkListing(fileName)
%
% input:
%   fileName:   path to the .egslst file of an EDKnrc run with conical
%               geometry (IFULL = ENERGY DEPOSITION KERNEL)
%
% output:
%   edk:        struct with fields
%                 energy          incident photon energy [MeV]
%                 numHistories    number of histories
%                 angleEdges      polar cone edges [deg], 1 x (nCones+1)
%                 radiusEdges     radial shell edges [mm], 1 x (nShells+1)
%                 total           fraction of incident energy deposited
%                                 per cone and shell, nCones x nShells
%                 primary         primary component, nCones x nShells
%                 scatter         total - primary, nCones x nShells
%                 totalRelUnc     relative 1 sigma uncertainty of total
%                 primaryRelUnc   relative 1 sigma uncertainty of primary
%                 scatterRelUnc   upper bound of the relative 1 sigma
%                                 uncertainty of scatter
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
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

txt = fileread(fileName);
lines = regexp(txt, '\r?\n', 'split');

edk = struct();
edk.energy = matRad_readScalar(txt, 'KINETIC ENERGY OF THE INCIDENT BEAM:\s*([\d.Ee+-]+)\s*MeV', NaN);
edk.numHistories = matRad_readScalar(txt, '(\d+)\s+histories run in', NaN);

headerIx = find(~cellfun(@isempty, strfind(lines, 'angle/deg.')), 1, 'first');
if isempty(headerIx)
    matRad_cfg.dispError('No conical energy deposition kernel table found in %s.', fileName);
end

rowPattern = '^\s*([\d.]+)\s+([\d.]+)\s+([\d.Ee+-]+)\s*\[\s*([\d.]+)%\s*\]\s+([\d.Ee+-]+)\s*\[\s*([\d.]+)%\s*\]';
data = zeros(numel(lines), 6);
numRows = 0;
for i = headerIx + 2:numel(lines)
    tok = regexp(lines{i}, rowPattern, 'tokens', 'once');
    if isempty(tok)
        break
    end
    numRows = numRows + 1;
    data(numRows, :) = str2double(tok);
end
data = data(1:numRows, :);

angles = unique(data(:, 1))';
radii = unique(data(:, 2))';
nCones = numel(angles);
nShells = numel(radii);

if nCones * nShells ~= numRows
    matRad_cfg.dispError('Unexpected table size in %s: %d rows for %d cones and %d shells.', fileName, numRows, nCones, nShells);
end

edk.angleEdges = [0 angles];
edk.radiusEdges = 10 * [0 radii];

[~, ic] = ismember(data(:, 1), angles);
[~, ir] = ismember(data(:, 2), radii);
ix = sub2ind([nCones nShells], ic, ir);

edk.total = zeros(nCones, nShells);
edk.primary = zeros(nCones, nShells);
edk.totalRelUnc = zeros(nCones, nShells);
edk.primaryRelUnc = zeros(nCones, nShells);

edk.total(ix) = data(:, 3);
edk.totalRelUnc(ix) = data(:, 4) / 100;
edk.primary(ix) = data(:, 5);
edk.primaryRelUnc(ix) = data(:, 6) / 100;

edk.scatter = max(edk.total - edk.primary, 0);
absUnc = sqrt((edk.totalRelUnc .* edk.total).^2 + (edk.primaryRelUnc .* edk.primary).^2);
edk.scatterRelUnc = inf(nCones, nShells);
hasScatter = edk.scatter > 0;
edk.scatterRelUnc(hasScatter) = absUnc(hasScatter) ./ edk.scatter(hasScatter);

end

function val = matRad_readScalar(txt, pattern, default)
tok = regexp(txt, pattern, 'tokens', 'once');
if isempty(tok)
    val = default;
else
    val = str2double(tok{1});
end
end
