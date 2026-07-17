function [cst, ixRing] = matRad_createRing(ct, cst, voiBase, vOuterMargin, vInnerMargin, varargin)
% matRad function to create a ring VOI around a base VOI, clipped to a
% limiting VOI
%
% The ring is obtained by growing the base VOI by the outer and the inner
% margin and subtracting the two. Growth is performed within the limiting
% VOI only, i.e. the ring cannot spread through structures outside of it.
% By default the ring is limited by all VOIs of the cst, i.e. by the patient
% outline. Note that the shape of the margin follows the connectivity used by
% matRad_addMargin and is therefore not a perfect box or sphere.
%
% By default each CT scenario is treated separately. With 'unionScenarios'
% the ring is instead grown around the union of the base VOI over all
% scenarios, e.g. to obtain a ring around the full motion of a 4D target.
% The limiting VOI is always applied per scenario, so the resulting ring may
% still differ between scenarios.
%
% All optional arguments can be given either as Name/Value pairs or as a
% metadata struct whose fields carry the same names, or as a combination of
% both (later values win). The struct thus mirrors the fifth column of the
% cst, so cst{i,5} of an existing VOI can be handed over directly.
%
% call:
%   [cst, ixRing] = matRad_createRing(ct, cst, voiBase, vOuterMargin, vInnerMargin)
%   [cst, ixRing] = matRad_createRing(ct, cst, voiBase, vOuterMargin, vInnerMargin, metadata)
%   [cst, ixRing] = matRad_createRing(ct, cst, voiBase, vOuterMargin, vInnerMargin, Name, Value)
%   [cst, ixRing] = matRad_createRing(ct, cst, voiBase, vOuterMargin, vInnerMargin, metadata, Name, Value)
%
% input:
%   ct:             matRad ct struct
%   cst:            matRad cst struct
%   voiBase:        base VOI the ring is grown around, either its name or its
%                   row index in the cst
%   vOuterMargin:   outer margin in mm, with fields x, y and z
%   vInnerMargin:   inner margin in mm, with fields x, y and z
%
%   Optional Name/Value properties (or fields of a metadata struct):
%   voiLimit:               VOI the ring is clipped to, either its name or its
%                           row index in the cst.
%                           Default: all VOIs of the cst, i.e. the patient outline
%   name:                   name of the ring VOI (cst column 2)
%                           Default: [name of the base VOI '_RING']
%   type:                   type of the ring VOI (cst column 3). Default: 'OAR'
%   objectives:             objective or cell array of objectives for the ring
%                           VOI (cst column 6). Default: {}
%   visibleColor:           color of the ring VOI. Default: [0 1 0]
%   Priority, Visible, TissueClass, alphaX, betaX:
%                           cst column 5 properties of the ring VOI. Each
%                           defaults to the value of the base VOI
%   diagonalConnectivity:   if true 26-connectivity is used for the margin,
%                           otherwise 6-connectivity. Default: false
%   unionScenarios:         if true the ring is grown around the union of the
%                           base VOI over all CT scenarios. Default: false
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

if nargin < 5
    matRadCfg.dispError('Not enough input arguments specified for matRad_createRing.');
end

% the base VOI is resolved first, the default ring name derives from it
ixBase = resolveVoiIx(cst, voiBase, 'voiBase', matRadCfg);

% StructExpand (default true) lets a metadata struct stand in for Name/Value
% pairs, so cst{i,5}-like structs and Name/Value pairs are interchangeable.
% PartialMatching is set explicitly because its default differs between
% Matlab (true) and Octave (false)
p = inputParser();
p.FunctionName = 'matRad_createRing';
p.StructExpand = true;
p.PartialMatching = false;
p.addParameter('voiLimit', []);
p.addParameter('name', [cst{ixBase, 2} '_RING']);
p.addParameter('type', 'OAR');
p.addParameter('objectives', {});
p.addParameter('visibleColor', [0 1 0]);
p.addParameter('Priority', []);
p.addParameter('Visible', []);
p.addParameter('TissueClass', []);
p.addParameter('alphaX', []);
p.addParameter('betaX', []);
p.addParameter('diagonalConnectivity', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
p.addParameter('unionScenarios', false, @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
p.parse(varargin{:});

metadata   = p.Results;
bDiaElem   = logical(metadata.diagonalConnectivity);
bUnionScen = logical(metadata.unionScenarios);

% without an explicit limiting VOI the ring is clipped to the patient outline,
% i.e. to the union of all VOIs of the cst
if isempty(metadata.voiLimit)
    limitRows = 1:size(cst, 1);
else
    limitRows = resolveVoiIx(cst, metadata.voiLimit, 'voiLimit', matRadCfg);
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

if numel(cst{ixBase, 4}) < ct.numOfCtScen
    matRadCfg.dispError('The base VOI ''%s'' must contain all CT scenarios.', cst{ixBase, 2});
end

if any(cellfun(@numel, cst(limitRows, 4)) < ct.numOfCtScen)
    matRadCfg.dispError('The limiting VOI(s) must contain all CT scenarios.');
end

% margins are grown by matRad_addMargin, which silently ignores negative
% values, so reject them here instead of returning a surprising ring
vOuter = [vOuterMargin.x vOuterMargin.y vOuterMargin.z];
vInner = [vInnerMargin.x vInnerMargin.y vInnerMargin.z];

if any(vOuter < 0) || any(vInner < 0)
    matRadCfg.dispError('Ring margins must not be negative.');
end

if any(vInner > vOuter)
    matRadCfg.dispError('The inner ring margin must not exceed the outer ring margin in any dimension.');
end

if any(strcmp(cst(:, 2), metadata.name))
    matRadCfg.dispError('A VOI named ''%s'' already exists in the cst.', metadata.name);
end

voiRing = cell(1, ct.numOfCtScen);

if bUnionScen
    % force column vectors, the cst does not guarantee an orientation
    baseScen = cellfun(@(v) v(:), cst{ixBase, 4}(1:ct.numOfCtScen), 'UniformOutput', false);
    baseUnion = unique(vertcat(baseScen{:}));
    matRadCfg.dispInfo('Ring VOI ''%s'' is grown around the union of ''%s'' over %d CT scenarios.\n', metadata.name, cst{ixBase, 2}, ct.numOfCtScen);
end

for scenIx = 1:ct.numOfCtScen
    if bUnionScen
        baseVoxels = baseUnion;
    else
        baseVoxels = cst{ixBase, 4}{scenIx};
    end

    % a single limiting VOI or, by default, the union of all VOIs of this scenario
    limitScen = cellfun(@(c) c{scenIx}(:), cst(limitRows, 4), 'UniformOutput', false);
    limitVoxels = unique(vertcat(limitScen{:}));

    baseMask = false(ct.cubeDim);
    baseMask(baseVoxels) = true;

    % matRad_addMargin grows the mask geodesically within the union of all
    % VOIs of the cst it is handed. Pass only the limiting VOI of the
    % current scenario so that the margin is confined to it right away.
    limitCst = cell(1, 6);
    limitCst{1, 4} = {limitVoxels};

    outerVoxels = find(matRad_addMargin(baseMask, limitCst, ct.resolution, vOuterMargin, bDiaElem));
    innerVoxels = find(matRad_addMargin(baseMask, limitCst, ct.resolution, vInnerMargin, bDiaElem));

    % the base VOI itself seeds the growth unconditionally and may reach
    % outside the limiting VOI, so clip the result once more
    voiRing{scenIx} = intersect(setdiff(outerVoxels, innerVoxels), limitVoxels);
end

ixRing = size(cst, 1) + 1;
ringId = max([cst{:, 1}]) + 1;

cst{ixRing, 1} = ringId;
cst{ixRing, 2} = metadata.name;
cst{ixRing, 3} = metadata.type;
cst{ixRing, 4} = voiRing;
% the ring starts from the cst column 5 properties of the base VOI, an
% explicitly requested property overrides the inherited one
cst{ixRing, 5} = cst{ixBase, 5};
cst{ixRing, 5}.visibleColor = metadata.visibleColor;

inheritedProps = {'Priority', 'Visible', 'TissueClass', 'alphaX', 'betaX'};
for propIx = 1:numel(inheritedProps)
    if ~isempty(metadata.(inheritedProps{propIx}))
        cst{ixRing, 5}.(inheritedProps{propIx}) = metadata.(inheritedProps{propIx});
    end
end

if isempty(metadata.Priority)
    if isfield(cst{ixBase, 5}, 'Priority')
        matRadCfg.dispInfo('Ring VOI ''%s'' inherits overlap priority %g from ''%s''.\n', metadata.name, cst{ixBase, 5}.Priority, cst{ixBase, 2});
    else
        matRadCfg.dispWarning('Neither the ring metadata nor the base VOI ''%s'' define an overlap priority for ''%s''.', cst{ixBase, 2}, metadata.name);
    end
end

if iscell(metadata.objectives)
    cst{ixRing, 6} = metadata.objectives;
else
    cst{ixRing, 6} = {metadata.objectives};
end

end

function ix = resolveVoiIx(cst, voi, argName, matRadCfg)
% resolves a VOI given as a cst row index or as a VOI name to a cst row index

% Matlab string scalars are accepted as names, on Octave isa simply returns false
if isa(voi, 'string') && isscalar(voi)
    voi = char(voi);
end

if isnumeric(voi) && isscalar(voi) && mod(voi, 1) == 0
    ix = voi;
    if ix < 1 || ix > size(cst, 1)
        matRadCfg.dispError('%s: index %d does not refer to an existing cst row (1 to %d).', argName, ix, size(cst, 1));
    end
    return;
end

if ~ischar(voi)
    matRadCfg.dispError('%s must be a VOI name or a cst row index.', argName);
end

ix = find(strcmp(cst(:, 2), voi));

if isempty(ix)
    matRadCfg.dispError('%s: the cst contains no VOI named ''%s''.', argName, voi);
elseif ~isscalar(ix)
    matRadCfg.dispError('%s: the cst contains %d VOIs named ''%s''.', argName, numel(ix), voi);
end

end
