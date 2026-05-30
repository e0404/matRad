function [cst, ixRing] = matRad_createRing(ixBase, ixLimit, cst, ct, vOuterMargin, vInnerMargin, metadata)
% matRad function to create an isotropic ring VOI clipped to a limiting VOI
%
% call:
%   [cst, ixRing] = matRad_createRing(ixBase, ixLimit, cst, ct, vOuterMargin, vInnerMargin, metadata)
%
% input:
%   ixBase:         row index of the base VOI in the cst struct
%   ixLimit:        row index of the limiting VOI in the cst struct
%   cst:            matRad cst struct
%   ct:             matRad ct struct
%   vOuterMargin:   outer margin in mm, with fields x, y and z
%   vInnerMargin:   inner margin in mm, with fields x, y and z
%   metadata:       struct with fields name, type and visibleColor for the
%                   created ring VOI
%
% output:
%   cst:            updated matRad cst struct
%   ixRing:         row index of the created ring VOI
%
% References
%   -
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

matRadCfg = MatRad_Config.instance();

if nargin < 7
    matRadCfg.dispError('Not enough input arguments specified for matRad_createRing.');
end

if ixBase < 1 || ixBase > size(cst, 1) || ixLimit < 1 || ixLimit > size(cst, 1)
    matRadCfg.dispError('Base and limiting VOI indices must refer to existing cst rows.');
end

if ~isfield(ct, 'numOfCtScen') || ct.numOfCtScen < 1
    matRadCfg.dispError('ct.numOfCtScen must be available and positive.');
end

if ~isfield(ct, 'cubeDim') || ~isfield(ct, 'resolution')
    matRadCfg.dispError('ct.cubeDim and ct.resolution are required to create a ring VOI.');
end

if ~all(isfield(vOuterMargin, {'x', 'y', 'z'})) || ~all(isfield(vInnerMargin, {'x', 'y', 'z'}))
    matRadCfg.dispError('Ring margins must contain x, y and z fields.');
end

if ~all(isfield(metadata, {'name', 'type', 'visibleColor'}))
    matRadCfg.dispError('Ring metadata must contain name, type and visibleColor fields.');
end

if numel(cst{ixBase, 4}) < ct.numOfCtScen || numel(cst{ixLimit, 4}) < ct.numOfCtScen
    matRadCfg.dispError('Base and limiting VOIs must contain all CT scenarios.');
end

voiRing = cell(1, ct.numOfCtScen);
useDiagonalConnectivity = false;

for scenIx = 1:ct.numOfCtScen
    baseVoxels = cst{ixBase, 4}{scenIx};
    limitVoxels = cst{ixLimit, 4}{scenIx};

    baseMask = zeros(ct.cubeDim);
    baseMask(baseVoxels) = 1;

    outerVoxels = find(matRad_addMargin(baseMask, cst, ct.resolution, vOuterMargin, useDiagonalConnectivity));
    innerVoxels = find(matRad_addMargin(baseMask, cst, ct.resolution, vInnerMargin, useDiagonalConnectivity));

    voiRing{scenIx} = intersect(setdiff(outerVoxels, innerVoxels), limitVoxels);
end

ixRing = size(cst, 1) + 1;

cst{ixRing, 1} = cst{end, 1} + 1;
cst{ixRing, 2} = metadata.name;
cst{ixRing, 3} = metadata.type;
cst{ixRing, 4} = voiRing;
cst{ixRing, 5} = cst{ixBase, 5};
cst{ixRing, 5}.visibleColor = metadata.visibleColor;

end
