function matRad_visApertureInfo(apertureInfo, mode, varargin)
% matRad function to visualize aperture shapes stored as struct
%
% call:
%   matRad_visApertureInfo(apertureInfo)
%   matRad_visApertureInfo(apertureInfo,mode)
%   matRad_visApertureInfo(apertureInfo,mode,Name,Value,...)
%
% input:
%   apertureInfo: aperture weight and shape info struct
%   mode:         view to show, either a single name or a cell array of
%                 names. Default: 'auto'
%                 'auto'       - 'perBeam' for static (IMRT/DAO) plans,
%                                'grid' + 'trajectory' + 'metrics' for VMAT
%                 'perBeam'    - one figure with one tab per beam, one
%                                subplot per shape of that beam. Suited to
%                                static delivery, where a beam carries many
%                                segments
%                 'grid'       - one (paginated) figure, one subplot per
%                                control point. Suited to VMAT, where every
%                                control point carries a single aperture
%                 'animate'    - single-aperture player, started playing:
%                                steps through all apertures in delivery
%                                order, each shown for a time proportional
%                                to its delivery time (VMAT) or weight. For
%                                VMAT the leaves and the gantry angle are
%                                interpolated between control points, so
%                                the aperture moves the way it is delivered.
%                                Playback holds the command line until it
%                                reaches the end of the plan, is paused, or
%                                the figure is closed
%                 'interactive'- the same player, started paused, so that
%                                the control point / shape sliders can be
%                                scrubbed by hand. Pausing an 'animate'
%                                figure leaves exactly this behind
%                 'trajectory' - leaf position vs gantry angle (VMAT only)
%                 'metrics'    - MU rate, gantry rotation speed and leaf
%                                speed vs gantry angle, against the machine
%                                delivery constraints (VMAT only)
%
%   Name-value pairs:
%   'leafCoordinate': how to place the leaf pairs on the vertical axis of
%                 the aperture plots, 'leafNum' (default) to label them by
%                 leaf pair index, or 'physical' to draw them at their
%                 physical coordinates
%
% output:
%   -
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

if nargin < 2 || isempty(mode)
    mode = 'auto';
end

[mode, varargin] = matRad_upgradeVisCall(mode, varargin);
opt = matRad_parseVisOptions(varargin);
leafCoord = opt.leafCoordinate;

% accept aperture info saved by an older matRad version
apertureInfo = matRad_upgradeApertureInfo(apertureInfo);

isVMAT = isfield(apertureInfo, 'runVMAT') && apertureInfo.runVMAT;

% resolve 'auto': the per-beam layout tiles the shapes of a single beam,
% which is the right thing for static delivery (few beams, many segments
% each) but degenerates for VMAT (many control points, one aperture each -
% one tab per control point). VMAT therefore gets the control-point grid
% plus the two views that actually show the arc: leaf trajectories and the
% delivery metrics against the machine constraints.
if iscell(mode)
    viewList = mode;
elseif strcmp(mode, 'auto')
    if isVMAT
        viewList = {'grid', 'trajectory', 'metrics'};
    else
        viewList = {'perBeam'};
    end
else
    viewList = {mode};
end

% global parameters
numOfBeams = size(apertureInfo.beam, 2);
wMax = matRad_getMaxWeight(apertureInfo, isVMAT);

for v = 1:numel(viewList)
    switch viewList{v}
        case 'perBeam'
            matRad_visPerBeam(apertureInfo, leafCoord, wMax, isVMAT, numOfBeams);
        case 'grid'
            matRad_visGrid(apertureInfo, leafCoord, wMax, numOfBeams);
        case {'animate', 'interactive'}
            matRad_visAperturePlayer(apertureInfo, leafCoord, wMax, isVMAT, numOfBeams, ...
                                     strcmp(viewList{v}, 'animate'));
        case 'trajectory'
            if matRad_requireVMAT(isVMAT, 'trajectory')
                matRad_visTrajectory(apertureInfo, numOfBeams);
            end
        case 'metrics'
            if matRad_requireVMAT(isVMAT, 'metrics')
                matRad_visMetrics(apertureInfo, numOfBeams);
            end
        otherwise
            matRad_cfg.dispError('Unknown view ''%s'' for matRad_visApertureInfo!', viewList{v});
    end
end

end

function [mode, args] = matRad_upgradeVisCall(mode, args)
% Accepts the old signature matRad_visApertureInfo(apertureInfo, leafCoord,
% viewType), where the second argument selected the leaf coordinates rather
% than the view. Those two values are not valid view names, so the old call
% can be told apart from the current one without ambiguity.

leafCoordinates = {'leafNum', 'physical'};

if ~ischar(mode) || ~any(strcmpi(mode, leafCoordinates))
    return
end

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispDeprecationWarning(['matRad_visApertureInfo(apertureInfo,''%s'',viewType) is ' ...
                                   'deprecated, use matRad_visApertureInfo(apertureInfo,view,' ...
                                   '''leafCoordinate'',''%s'') instead.'], mode, mode);

legacyLeafCoord = mode;

% the old third argument, if any, was the view
mode = 'auto';
if ~isempty(args) && (ischar(args{1}) || iscell(args{1}))
    mode = args{1};
    args = args(2:end);
end

args = [{'leafCoordinate', legacyLeafCoord}, args];
end

function opt = matRad_parseVisOptions(args)
% Name-value options of the aperture visualization.

matRad_cfg = MatRad_Config.instance();

opt.leafCoordinate = 'leafNum';
leafCoordinates = {'leafNum', 'physical'};

if mod(numel(args), 2) ~= 0
    matRad_cfg.dispError('The optional arguments of matRad_visApertureInfo have to be name-value pairs!');
end

for k = 1:2:numel(args)

    name = args{k};
    value = args{k + 1};

    if ~ischar(name)
        matRad_cfg.dispError('Option names of matRad_visApertureInfo have to be character arrays!');
    end

    switch lower(name)
        case 'leafcoordinate'
            ix = find(strcmpi(value, leafCoordinates), 1);
            if isempty(ix)
                matRad_cfg.dispError('leafCoordinate has to be one of ''%s''!', ...
                                     strjoin(leafCoordinates, ''', '''));
            end
            opt.leafCoordinate = leafCoordinates{ix};
        otherwise
            matRad_cfg.dispError('Unknown option ''%s'' for matRad_visApertureInfo!', name);
    end
end
end

function ok = matRad_requireVMAT(isVMAT, viewName)
% The trajectory and metrics views read apertureInfo.arc, which only
% exists for VMAT plans.
ok = isVMAT;
if ~ok
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispWarning('The ''%s'' view is only available for VMAT plans, skipping it.', viewName);
end
end

function wMax = matRad_getMaxWeight(apertureInfo, isVMAT)
% For VMAT the weight colouring is normalized across the whole arc, so that
% shading is comparable between control points. For static delivery it is
% normalized per beam (handled in matRad_visPerBeam) - here we only provide the
% global fallback.
wMax = 0;
for i = 1:size(apertureInfo.beam, 2)
    for j = 1:numel(apertureInfo.beam(i).shape)
        wMax = max(wMax, apertureInfo.beam(i).shape(j).weight);
    end
end
if wMax <= 0 || ~isVMAT && wMax == 0
    wMax = 1; % degenerate plan, avoid division by zero
end
end

%% shared helpers of the aperture views

function n = matRad_apertureNumShapes(apertureInfo, i)
% Deliberately count the shape structs rather than reading
% beam(i).numOfShapes: in VMAT, interpolated (non-DAO) beams carry a
% computed shape while their numOfShapes (the number of *optimized* shapes)
% is 0 - and those interpolated apertures should be shown too.
n = numel(apertureInfo.beam(i).shape);
end

function panels = matRad_aperturePanels(apertureInfo, numOfBeams)
% Flattens all (beam, shape) pairs into one list, in delivery order: the
% beams are kept in the order they are stored in, which for VMAT is the
% order in which the arc is traversed (sorting by gantry angle would break
% reverse arcs and multi-arc plans).

panels = zeros(0, 2);
for i = 1:numOfBeams
    for j = 1:matRad_apertureNumShapes(apertureInfo, i)
        panels(end + 1, :) = [i j]; %#ok<AGROW>
    end
end
end

function isDAOBeam = matRad_apertureDAOFlags(apertureInfo, numOfBeams)
% Only VMAT distinguishes DAO control points from interpolated ones; for any
% other plan every beam counts as optimized.
isDAOBeam = true(1, numOfBeams);
if isfield(apertureInfo, 'arc') && isfield(apertureInfo.arc, 'beam')
    isDAOBeam = logical([apertureInfo.arc.beam.isDAOBeam]);
end
end

function val = matRad_apertureFieldOrNaN(s, fieldName)
% Which fields an aperture info struct carries depends on how far a plan has
% been taken - the gantry angle and the delivery timing, for instance, are
% VMAT-only. Reading them through this lets a caller display what is there
% and leave out what is not.
val = NaN;
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName)(1);
end
end

function c = matRad_apertureWeightColor(weight, wMax)
% Background colour encoding an aperture weight: light grey for an
% unweighted aperture, saturated red for the heaviest one.
color = [0.2:0.01:0.8; 0.2:0.01:0.8; 0.2:0.01:0.8]';
color = flipud(color);
color(:, 3) = 0;
color(:, 2) = 0;

colorInd = max(ceil((weight / wMax) * (size(color, 1) - 1) + eps), 1);
colorInd = min(colorInd, size(color, 1));
c = color(colorInd, :);
end

function [rowLo, rowHi] = matRad_leafRowBounds(apertureInfo, i, leafCoord, bixelWidth)
% Lower/upper edge of every active leaf row of beam i, as row vectors.

nK = apertureInfo.beam(i).numOfActiveLeafPairs;

if strcmp(leafCoord, 'physical')
    rowCentre = apertureInfo.beam(i).leafPairPos(1:nK);
    halfRow = bixelWidth / 2;
else
    % the leaf indices have to be flipped in order to fit to the order of
    % the leaf positions (1st row of leafPos is lowest row in physical
    % coordinates)
    rowCentre = flipud(find(apertureInfo.beam(i).isActiveLeafPair));
    rowCentre = rowCentre(1:nK);
    halfRow = 1 / 2;
end

rowCentre = reshape(rowCentre, 1, []);
rowLo = rowCentre - halfRow;
rowHi = rowCentre + halfRow;
end

function [xLeft, xRight, yRow] = matRad_apertureQuads(apertureInfo, i, leftLeafPos, rightLeafPos, leafCoord, bixelWidth)
% Vertices of the closed (blocked) parts of the leaf pairs, as 4-by-nLeafPairs
% matrices - one column per quadrangle, i.e. one patch/fill face per leaf.
% Taking the leaf positions as an argument rather than reading them from a
% shape lets the player draw positions interpolated between two shapes.

minX = apertureInfo.beam(i).MLCWindow(1);
maxX = apertureInfo.beam(i).MLCWindow(2);
nK = apertureInfo.beam(i).numOfActiveLeafPairs;

[rowLo, rowHi] = matRad_leafRowBounds(apertureInfo, i, leafCoord, bixelWidth);
yRow = [rowLo; rowLo; rowHi; rowHi];

leftLeafPos = reshape(leftLeafPos(1:nK), 1, []);
rightLeafPos = reshape(rightLeafPos(1:nK), 1, []);

xLeft = [repmat(minX, 1, nK); leftLeafPos; leftLeafPos; repmat(minX, 1, nK)];
xRight = [rightLeafPos; repmat(maxX, 1, nK); repmat(maxX, 1, nK); rightLeafPos];
end

function matRad_labelApertureAxes(hAx, leafCoord)
xlabel(hAx, 'horiz. pos. [mm]');
if strcmp(leafCoord, 'physical')
    ylabel(hAx, 'vert. pos. [mm]');
else
    ylabel(hAx, 'leaf pair #');
end
end

function matRad_drawAperture(hAx, apertureInfo, i, j, leafCoord, bixelWidth)
% Draws the closed leaf pairs of shape j of beam i into the given axes.

[xLeft, xRight, yRow] = matRad_apertureQuads(apertureInfo, i, ...
                                             apertureInfo.beam(i).shape(j).leftLeafPos, ...
                                             apertureInfo.beam(i).shape(j).rightLeafPos, ...
                                             leafCoord, bixelWidth);

% left bank, then right bank - one fill object each, one face per leaf pair
fill(hAx, xLeft, yRow, [0.5 0.5 0.5]);
fill(hAx, xRight, yRow, [0.5 0.5 0.5]);

axis(hAx, 'tight');
matRad_labelApertureAxes(hAx, leafCoord);
end

function matRad_visPerBeam(apertureInfo, leafCoord, wMaxGlobal, isVMAT, numOfBeams)
% One tab per beam, one subplot per shape of that beam. Where tabs are not
% available (Octave does not implement uitabgroup) this falls back to the
% original layout of one figure per beam.

bixelWidth = apertureInfo.bixelWidth;

% the shared figure has to fit the busiest beam
maxShapes = 0;
for i = 1:numOfBeams
    maxShapes = max(maxShapes, matRad_apertureNumShapes(apertureInfo, i));
end
[maxLines, maxColumns] = matRad_perBeamLayout(maxShapes);

hTabGroup = matRad_tryTabGroup([1.8 * maxColumns 3 * maxLines], 'Aperture shapes');

for i = 1:numOfBeams

    nShapes = matRad_apertureNumShapes(apertureInfo, i);

    if isVMAT
        wMax = wMaxGlobal;
    else
        % if not VMAT, let wMax be the max weight of a particular angle
        wMax = wMaxGlobal;
        if nShapes > 0
            wMax = max([apertureInfo.beam(i).shape(:).weight]);
        end
    end

    [subplotLines, subplotColumns] = matRad_perBeamLayout(nShapes);

    if isempty(hTabGroup)
        % open a new figure for every beam
        hParent = figure('units', 'inches');
        set(hParent, 'pos', [0 0 1.8 * subplotColumns 3 * subplotLines], ...
            'Name', matRad_beamLabel(apertureInfo, i), 'NumberTitle', 'off');
    else
        hParent = uitab('Parent', hTabGroup, 'Title', matRad_beamLabel(apertureInfo, i));
    end

    % loop over all shapes of the beam
    for j = 1:nShapes

        hAx = subplot(subplotLines, subplotColumns, j, 'Parent', hParent);

        title(hAx, ['Beam: ' num2str(i) ' Shape: ' num2str(j) ' w=' ...
                    num2str(apertureInfo.beam(i).shape(j).weight, 2)], ...
              'Fontsize', 8);

        set(hAx, 'Color', matRad_apertureWeightColor(apertureInfo.beam(i).shape(j).weight, wMax));
        hold(hAx, 'on');

        matRad_drawAperture(hAx, apertureInfo, i, j, leafCoord, bixelWidth);
    end
end
end

function [subplotLines, subplotColumns] = matRad_perBeamLayout(nShapes)
% Two rows of shapes, as wide as it takes.
subplotColumns = max(ceil(nShapes / 2), 1);
subplotLines = max(ceil(nShapes / subplotColumns), 1);
end

function label = matRad_beamLabel(apertureInfo, i)
angle = matRad_apertureFieldOrNaN(apertureInfo.beam(i), 'gantryAngle');
if isfinite(angle)
    label = sprintf('Beam %d (%.0f%s)', i, angle, char(176));
else
    label = sprintf('Beam %d', i);
end
end

function hTabGroup = matRad_tryTabGroup(figSizeInches, figName)
% Creates a figure holding a tab group, or returns empty where tabs are not
% supported, leaving it to the caller to fall back to separate figures.

hFig = [];

try
    hFig = figure('units', 'inches', 'Name', figName, 'NumberTitle', 'off');
    set(hFig, 'pos', [0 0 max(figSizeInches(1), 4) max(figSizeInches(2), 3)]);
    hTabGroup = uitabgroup('Parent', hFig);
catch
    if ~isempty(hFig) && ishandle(hFig)
        close(hFig);
    end
    hTabGroup = [];
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispDebug(['uitabgroup is not available in this environment, ' ...
                          'showing one figure per beam instead.\n']);
end
end

function matRad_visGrid(apertureInfo, leafCoord, wMax, numOfBeams)
% One subplot per control point, paginated over as many figures as needed.
% This is the VMAT counterpart of matRad_visPerBeam: there the shapes of one beam
% are tiled, here the (single) aperture of every beam is.

bixelWidth = apertureInfo.bixelWidth;
maxPerFigure = 30; % a 6x5 grid is still legible, and fits a typical single arc

panels = matRad_aperturePanels(apertureInfo, numOfBeams);
numPanels = size(panels, 1);
if numPanels == 0
    return
end

% Spread the panels evenly over the pages instead of filling each to the
% brim, so that e.g. 31 apertures give two balanced pages rather than a
% full one followed by a page holding a single aperture.
numPages = ceil(numPanels / maxPerFigure);
perPage = ceil(numPanels / numPages);
isDAOBeam = matRad_apertureDAOFlags(apertureInfo, numOfBeams);

for page = 1:numPages

    first = (page - 1) * perPage + 1;
    last = min(page * perPage, numPanels);
    nThisPage = last - first + 1;

    subplotColumns = ceil(sqrt(nThisPage));
    subplotLines = ceil(nThisPage / subplotColumns);

    figure('units', 'inches');
    set(gcf, 'pos', [0 0 1.8 * subplotColumns 2.2 * subplotLines]);
    if numPages > 1
        set(gcf, 'Name', sprintf('Apertures %d-%d of %d', first, last, numPanels));
    else
        set(gcf, 'Name', 'Apertures');
    end

    for p = first:last

        i = panels(p, 1);
        j = panels(p, 2);

        hAx = subplot(subplotLines, subplotColumns, p - first + 1);

        % Label by gantry angle rather than beam index where it is known -
        % along an arc the angle is what identifies a control point.
        % Interpolated (non-DAO) control points are marked, since their
        % apertures are derived rather than optimized.
        angle = matRad_apertureFieldOrNaN(apertureInfo.beam(i), 'gantryAngle');
        if isfinite(angle)
            panelLabel = sprintf('G=%.0f%s', angle, char(176));
        else
            panelLabel = sprintf('Beam %d/%d', i, j);
        end
        if ~isDAOBeam(i)
            panelLabel = [panelLabel ' (interp)'];
        end
        title(hAx, sprintf('%s w=%s', panelLabel, ...
                           num2str(apertureInfo.beam(i).shape(j).weight, 2)), 'Fontsize', 7);

        set(hAx, 'Color', matRad_apertureWeightColor(apertureInfo.beam(i).shape(j).weight, wMax));
        hold(hAx, 'on');

        matRad_drawAperture(hAx, apertureInfo, i, j, leafCoord, bixelWidth);
        set(hAx, 'Fontsize', 6);
    end
end
end

function matRad_visTrajectory(apertureInfo, numOfBeams)
% Leaf position vs gantry angle - the standard VMAT sweep plot. Makes the
% leaf trajectory, the sweep direction reversals between arc sectors and
% any leaf-touching artifacts directly visible.

gantryAngles = [apertureInfo.beam.gantryAngle];
[gantryAngles, order] = sort(gantryAngles);

numLeafPairs = apertureInfo.beam(1).numOfActiveLeafPairs;
leftPos = nan(numLeafPairs, numOfBeams);
rightPos = nan(numLeafPairs, numOfBeams);

for b = 1:numOfBeams
    i = order(b);
    if matRad_apertureNumShapes(apertureInfo, i) < 1
        continue
    end
    leftPos(:, b) = apertureInfo.beam(i).shape(1).leftLeafPos;
    rightPos(:, b) = apertureInfo.beam(i).shape(1).rightLeafPos;
end

figure('units', 'inches');
set(gcf, 'pos', [0 0 9 5], 'Name', 'Leaf trajectories');
hold on;

% mark the DAO control points: leaves are optimized there and merely
% interpolated in between
isDAOBeam = matRad_apertureDAOFlags(apertureInfo, numOfBeams);
daoAngles = [apertureInfo.beam(isDAOBeam).gantryAngle];
yLo = min([leftPos(:); rightPos(:)]);
yHi = max([leftPos(:); rightPos(:)]);
for a = 1:numel(daoAngles)
    plot([daoAngles(a) daoAngles(a)], [yLo yHi], ':', 'Color', [0.8 0.8 0.8]);
end

hLeft = [];
hRight = [];
for k = 1:numLeafPairs
    hLeft = plot(gantryAngles, leftPos(k, :), '-', 'Color', [0.2 0.4 0.8], 'LineWidth', 0.5);
    hRight = plot(gantryAngles, rightPos(k, :), '-', 'Color', [0.8 0.3 0.2], 'LineWidth', 0.5);
end

% highlight the central leaf pair, which is usually the most instructive
kCentre = max(round(numLeafPairs / 2), 1);
hLeftC = plot(gantryAngles, leftPos(kCentre, :), '-', 'Color', [0.1 0.2 0.5], 'LineWidth', 2);
hRightC = plot(gantryAngles, rightPos(kCentre, :), '-', 'Color', [0.5 0.1 0.1], 'LineWidth', 2);

xlabel(['gantry angle [' char(176) ']']);
ylabel('leaf position [mm]');
title(sprintf('Leaf trajectories (%d leaf pairs, %d control points)', numLeafPairs, numOfBeams));
if ~isempty(hLeft) && ~isempty(hRight)
    legend([hLeft hRight hLeftC hRightC], ...
           {'left bank', 'right bank', ...
            sprintf('left, pair %d', kCentre), sprintf('right, pair %d', kCentre)}, ...
           'Location', 'best');
end
grid on;
box on;
end

function matRad_visMetrics(apertureInfo, numOfBeams)
% Delivery metrics along the arc, plotted against the machine constraints.
% Values live on the DAO control points only, since those are what the
% optimizer constrains.

isDAOBeam = matRad_apertureDAOFlags(apertureInfo, numOfBeams);

hasShape = false(1, numOfBeams);
for i = 1:numOfBeams
    hasShape(i) = matRad_apertureNumShapes(apertureInfo, i) >= 1;
end
beamIx = find(isDAOBeam & hasShape);

angles = zeros(1, numel(beamIx));
muRate = zeros(1, numel(beamIx));
gantryRot = zeros(1, numel(beamIx));
leafSpeed = zeros(1, numel(beamIx));

for b = 1:numel(beamIx)
    i = beamIx(b);
    angles(b) = apertureInfo.beam(i).gantryAngle;
    muRate(b) = matRad_apertureFieldOrNaN(apertureInfo.beam(i).shape(1), 'MURate');
    gantryRot(b) = matRad_apertureFieldOrNaN(apertureInfo.beam(i), 'gantryRot');
    leafSpeed(b) = matRad_apertureFieldOrNaN(apertureInfo.beam(i), 'maxLeafSpeed');
end

[angles, order] = sort(angles);
muRate = muRate(order);
gantryRot = gantryRot(order);
leafSpeed = leafSpeed(order);

constraints = struct();
if isfield(apertureInfo, 'deliveryConstraints')
    constraints = apertureInfo.deliveryConstraints;
end

figure('units', 'inches');
set(gcf, 'pos', [0 0 8 7], 'Name', 'VMAT delivery metrics');

matRad_metricPlot(1, angles, muRate, constraints, 'monitorUnitRate', 'MU rate [MU/s]');
matRad_metricPlot(2, angles, gantryRot, constraints, 'gantryRotationSpeed', ...
                  ['gantry speed [' char(176) '/s]']);
matRad_metricPlot(3, angles, leafSpeed, constraints, 'leafSpeed', 'max leaf speed [mm/s]');
end

function matRad_metricPlot(idx, angles, values, constraints, constraintName, label)
% One metric panel with its lower/upper machine limit drawn in.

subplot(3, 1, idx);
hold on;

plot(angles, values, 'o-', 'Color', [0.2 0.4 0.8], 'MarkerSize', 3);

if isfield(constraints, constraintName) && ~isempty(angles)
    lim = constraints.(constraintName);
    for b = 1:numel(lim)
        plot([min(angles) max(angles)], [lim(b) lim(b)], '--', 'Color', [0.8 0.2 0.2]);
    end
    % keep both the data and the limits visible
    yLo = min([values(:); lim(:)]);
    yHi = max([values(:); lim(:)]);
    if isfinite(yLo) && isfinite(yHi) && yHi > yLo
        ylim([yLo - 0.05 * (yHi - yLo), yHi + 0.05 * (yHi - yLo)]);
    end
end

xlabel(['gantry angle [' char(176) ']']);
ylabel(label);
grid on;
box on;
end

%% player (the 'animate' / 'interactive' views)

function matRad_visAperturePlayer(apertureInfo, leafCoord, wMax, isVMAT, numOfBeams, autoPlay)
% Single-aperture player: one axes showing one aperture at a time, driven
% either by the transport (play/pause, speed, loop) or by the two sliders
% (control point / beam, and shape within that beam). 'animate' and
% 'interactive' are the same figure, they only differ in whether it starts
% playing - pausing an animation leaves the sliders free to scrub.

matRad_cfg = MatRad_Config.instance();

data = struct();
data.apertureInfo = apertureInfo;
data.leafCoordinate = leafCoord;
data.isVMAT = isVMAT;
data.wMax = wMax;
data.bixelWidth = apertureInfo.bixelWidth;
data.numOfBeams = numOfBeams;
data.isDAOBeam = matRad_apertureDAOFlags(apertureInfo, numOfBeams);
data.panels = matRad_aperturePanels(apertureInfo, numOfBeams);

numPanels = size(data.panels, 1);
if numPanels == 0
    matRad_cfg.dispWarning('No apertures to show, skipping the aperture player.');
    return
end

[data.duration, data.timeIsPhysical] = matRad_panelDurations(apertureInfo, data.panels, isVMAT);
data.tStart = [0; cumsum(data.duration(1:end - 1))];
data.totalTime = sum(data.duration);

if data.timeIsPhysical
    % playback at 1x then runs in real delivery time
    data.timeLabel = 'delivery time';
    data.timeKind = 'delivery';
else
    % no machine timing available: the dwell of an aperture is taken
    % proportional to its weight, scaled so that the plan plays back in a
    % few seconds at 1x
    data.timeLabel = 'dwell';
    data.timeKind = 'playback';
end

% Fixed axis limits over the whole plan, so that the aperture does not
% jump around while the animation runs.
[data.xLim, data.yLim] = matRad_playerAxisLimits(apertureInfo, leafCoord, data.bixelWidth, numOfBeams);

hFig = matRad_playerFigure(data, isVMAT);

state = struct();
state.curTime = 0;
state.isPlaying = false;
state.speed = 1;
state.doLoop = false;
state.doInterp = isVMAT;
state.shownBeam = 0;
state.shownShape = 0;
setappdata(hFig, 'matRadApertureState', state);

matRad_playerDraw(hFig);

if autoPlay
    matRad_playerPlay(hFig);
end
end

function hFig = matRad_playerFigure(data, isVMAT)
% Builds the figure, its axes and all controls, and stores the static part
% of the player (aperture data and handles) in the figure's appdata.

hFig = figure('Units', 'pixels', 'Position', [100 100 760 700], ...
              'Name', 'Aperture player', 'NumberTitle', 'off', ...
              'Color', get(0, 'defaultUicontrolBackgroundColor'));

data.hAx = axes('Parent', hFig, 'Units', 'normalized', 'Position', [0.11 0.45 0.84 0.47]);
hold(data.hAx, 'on');
box(data.hAx, 'on');
set(data.hAx, 'XLim', data.xLim, 'YLim', data.yLim);
matRad_labelApertureAxes(data.hAx, data.leafCoordinate);
data.hTitle = title(data.hAx, '');

% one patch per leaf bank, holding one face per leaf pair - updating their
% vertices is what animates the aperture
i = data.panels(1, 1);
j = data.panels(1, 2);
[xLeft, xRight, yRow] = matRad_apertureQuads(data.apertureInfo, i, ...
                                             data.apertureInfo.beam(i).shape(j).leftLeafPos, ...
                                             data.apertureInfo.beam(i).shape(j).rightLeafPos, ...
                                             data.leafCoordinate, data.bixelWidth);
data.hLeft = patch('Parent', data.hAx, 'XData', xLeft, 'YData', yRow, ...
                   'FaceColor', [0.5 0.5 0.5]);
data.hRight = patch('Parent', data.hAx, 'XData', xRight, 'YData', yRow, ...
                    'FaceColor', [0.5 0.5 0.5]);

data.hInfo = matRad_playerText(hFig, [0.06 0.330 0.90 0.045], '', 'left');
set(data.hInfo, 'FontSize', 9);

% control point / beam slider
matRad_playerText(hFig, [0.03 0.250 0.18 0.04], 'control point', 'left');
data.hBeamSlider = matRad_playerSlider(hFig, [0.22 0.255 0.56 0.035], data.numOfBeams, ...
                                       {@matRad_playerSliderCb, 'beam'});
data.hBeamValue = matRad_playerText(hFig, [0.80 0.250 0.18 0.04], '', 'left');

% shape slider (VMAT carries a single shape per control point, so it stays
% disabled there)
matRad_playerText(hFig, [0.03 0.190 0.18 0.04], 'shape', 'left');
data.hShapeSlider = matRad_playerSlider(hFig, [0.22 0.195 0.56 0.035], ...
                                        matRad_maxNumShapes(data), ...
                                        {@matRad_playerSliderCb, 'shape'});
data.hShapeValue = matRad_playerText(hFig, [0.80 0.190 0.18 0.04], '', 'left');

% transport
data.hPlay = uicontrol(hFig, 'Style', 'pushbutton', 'Units', 'normalized', ...
                       'Position', [0.03 0.055 0.13 0.075], 'String', 'Play', ...
                       'Callback', @matRad_playerPlayCb);

matRad_playerText(hFig, [0.18 0.062 0.07 0.045], 'speed', 'left');
data.speedList = [0.25 0.5 1 2 4 8];
data.hSpeed = uicontrol(hFig, 'Style', 'popupmenu', 'Units', 'normalized', ...
                        'Position', [0.25 0.065 0.10 0.05], ...
                        'String', {'0.25x', '0.5x', '1x', '2x', '4x', '8x'}, ...
                        'Value', find(data.speedList == 1, 1), ...
                        'Callback', @matRad_playerSpeedCb);

% Looping is off by default: while the animation runs it holds the command
% line (it has to, to stay responsive to the buttons in Octave as well as
% MATLAB), so playing a plan once and returning is the friendlier default.
data.hLoop = uicontrol(hFig, 'Style', 'checkbox', 'Units', 'normalized', ...
                       'Position', [0.375 0.065 0.11 0.045], 'String', 'loop', ...
                       'Value', 0, 'Callback', @matRad_playerToggleCb);

data.hInterp = uicontrol(hFig, 'Style', 'checkbox', 'Units', 'normalized', ...
                         'Position', [0.49 0.065 0.24 0.045], 'String', 'dynamic leaves', ...
                         'Value', double(isVMAT), 'Enable', matRad_onOff(isVMAT), ...
                         'TooltipString', ['interpolate the leaf positions and the gantry ' ...
                                           'angle between consecutive control points'], ...
                         'Callback', @matRad_playerToggleCb);

data.hTime = matRad_playerText(hFig, [0.74 0.062 0.24 0.045], '', 'right');

setappdata(hFig, 'matRadApertureData', data);
end

function h = matRad_playerText(hFig, pos, str, align)
h = uicontrol(hFig, 'Style', 'text', 'Units', 'normalized', 'Position', pos, ...
              'String', str, 'HorizontalAlignment', align, ...
              'BackgroundColor', get(hFig, 'Color'));
end

function h = matRad_playerSlider(hFig, pos, numSteps, callback)
% A slider over 1:numSteps. Sliders need Min < Max even when there is
% nothing to slide over, so a degenerate range is created disabled.

h = uicontrol(hFig, 'Style', 'slider', 'Units', 'normalized', 'Position', pos, ...
              'Min', 1, 'Max', max(numSteps, 2), 'Value', 1, 'Callback', callback);

if numSteps > 1
    step = 1 / (numSteps - 1);
    set(h, 'SliderStep', [step, max(step, 0.1)]);
else
    set(h, 'Enable', 'off');
end

try
    % live update while dragging, rather than only on release
    addlistener(h, 'ContinuousValueChange', @(src, evt) callback{1}(src, evt, callback{2}));
catch
    % ContinuousValueChange is a MATLAB-only event - without it the slider
    % simply updates when it is released
end
end

function s = matRad_onOff(tf)
if tf
    s = 'on';
else
    s = 'off';
end
end

function n = matRad_maxNumShapes(data)
n = 1;
for i = 1:data.numOfBeams
    n = max(n, matRad_apertureNumShapes(data.apertureInfo, i));
end
end

function [duration, timeIsPhysical] = matRad_panelDurations(apertureInfo, panels, isVMAT)
% Dwell time of every aperture. VMAT carries the machine timing, so the
% animation can run in real delivery time; without it (static delivery, or
% an apertureInfo straight from the sequencer) the dwell is taken
% proportional to the shape weight instead.

numPanels = size(panels, 1);
duration = zeros(numPanels, 1);
timeIsPhysical = isVMAT;

for p = 1:numPanels
    i = panels(p, 1);
    j = panels(p, 2);

    t = NaN;
    if isVMAT
        t = matRad_apertureFieldOrNaN(apertureInfo.beam(i), 'time');
        if ~(isfinite(t) && t > 0)
            % before optimization only the gantry speed is set
            gantryRot = matRad_apertureFieldOrNaN(apertureInfo.beam(i), 'gantryRot');
            arcDiff = matRad_arcFieldOrNaN(apertureInfo, i, 'doseAngleBordersDiff');
            if isfinite(gantryRot) && gantryRot > 0 && isfinite(arcDiff)
                t = arcDiff / gantryRot;
            end
        end
    end

    if ~(isfinite(t) && t > 0)
        timeIsPhysical = false;
        t = apertureInfo.beam(i).shape(j).weight;
    end

    duration(p) = t;
end

duration(~isfinite(duration) | duration < 0) = 0;

if all(duration == 0)
    duration(:) = 1;
elseif any(duration == 0)
    % zero-weight apertures would otherwise be skipped entirely - give them
    % a short blip so that they are still visible in the sequence
    duration(duration == 0) = 0.05 * mean(duration(duration > 0));
end

if ~timeIsPhysical
    % playback of the whole plan takes ~0.4 s per aperture at 1x
    duration = duration * (0.4 * numPanels / sum(duration));
end
end

function val = matRad_arcFieldOrNaN(apertureInfo, i, fieldName)
val = NaN;
if isfield(apertureInfo, 'arc') && isfield(apertureInfo.arc, 'beam') && ...
        numel(apertureInfo.arc.beam) >= i
    val = matRad_apertureFieldOrNaN(apertureInfo.arc.beam(i), fieldName);
end
end

function [xLim, yLim] = matRad_playerAxisLimits(apertureInfo, leafCoord, bixelWidth, numOfBeams)

xLo = inf;
xHi = -inf;
yLo = inf;
yHi = -inf;

for i = 1:numOfBeams
    if matRad_apertureNumShapes(apertureInfo, i) < 1
        continue
    end
    xLo = min(xLo, apertureInfo.beam(i).MLCWindow(1));
    xHi = max(xHi, apertureInfo.beam(i).MLCWindow(2));

    [rowLo, rowHi] = matRad_leafRowBounds(apertureInfo, i, leafCoord, bixelWidth);
    yLo = min(yLo, min(rowLo));
    yHi = max(yHi, max(rowHi));
end

if ~(isfinite(xLo) && isfinite(xHi) && xHi > xLo)
    xLo = -1;
    xHi = 1;
end
if ~(isfinite(yLo) && isfinite(yHi) && yHi > yLo)
    yLo = 0;
    yHi = 1;
end

xLim = [xLo xHi];
yLim = [yLo yHi];
end

function [leftLeafPos, rightLeafPos, gantryAngle] = matRad_playerAperture(data, p, frac, doInterp)
% Aperture shown at fraction frac of the dwell of panel p. With dynamic
% leaves the aperture of a control point is taken as the machine state at
% the *start* of its sector, and leaves and gantry travel linearly from
% there to the next control point over the dwell - which is what the
% delivery constraints on leaf speed and gantry speed model as well.

i = data.panels(p, 1);
j = data.panels(p, 2);

leftLeafPos = data.apertureInfo.beam(i).shape(j).leftLeafPos(:);
rightLeafPos = data.apertureInfo.beam(i).shape(j).rightLeafPos(:);
gantryAngle = matRad_apertureFieldOrNaN(data.apertureInfo.beam(i), 'gantryAngle');

if ~doInterp || frac <= 0 || p >= size(data.panels, 1)
    return
end

iNext = data.panels(p + 1, 1);
jNext = data.panels(p + 1, 2);
leftNext = data.apertureInfo.beam(iNext).shape(jNext).leftLeafPos(:);
rightNext = data.apertureInfo.beam(iNext).shape(jNext).rightLeafPos(:);

% a change in the number of active leaf pairs leaves nothing to interpolate
if numel(leftNext) == numel(leftLeafPos) && numel(rightNext) == numel(rightLeafPos)
    leftLeafPos = (1 - frac) * leftLeafPos + frac * leftNext;
    rightLeafPos = (1 - frac) * rightLeafPos + frac * rightNext;
end

% take the short way round, so that e.g. 350 -> 10 degrees does not sweep
% backwards through the whole arc
dAngle = mod(matRad_apertureFieldOrNaN(data.apertureInfo.beam(iNext), 'gantryAngle') - gantryAngle + 180, 360) - 180;
gantryAngle = gantryAngle + frac * dAngle;
end

function [p, frac] = matRad_playerPanelAtTime(data, t)
p = find(data.tStart <= t, 1, 'last');
if isempty(p)
    p = 1;
end
frac = 0;
if data.duration(p) > 0
    frac = min(max((t - data.tStart(p)) / data.duration(p), 0), 1);
end
end

function matRad_playerDraw(hFig)
% Redraws the aperture for the current playback time and syncs all labels
% and sliders to it.

data = getappdata(hFig, 'matRadApertureData');
state = getappdata(hFig, 'matRadApertureState');

[p, frac] = matRad_playerPanelAtTime(data, state.curTime);
i = data.panels(p, 1);
j = data.panels(p, 2);

[leftLeafPos, rightLeafPos, gantryAngle] = matRad_playerAperture(data, p, frac, state.doInterp);
[xLeft, xRight, yRow] = matRad_apertureQuads(data.apertureInfo, i, leftLeafPos, rightLeafPos, ...
                                             data.leafCoordinate, data.bixelWidth);

set(data.hLeft, 'XData', xLeft, 'YData', yRow);
set(data.hRight, 'XData', xRight, 'YData', yRow);

weight = data.apertureInfo.beam(i).shape(j).weight;
set(data.hAx, 'Color', matRad_apertureWeightColor(weight, data.wMax));

nShapes = matRad_apertureNumShapes(data.apertureInfo, i);
titleStr = sprintf('Beam %d/%d, shape %d/%d', i, data.numOfBeams, j, nShapes);
if isfinite(gantryAngle)
    titleStr = sprintf('%s - G = %.1f%s', titleStr, gantryAngle, char(176));
end
set(data.hTitle, 'String', titleStr);
set(data.hInfo, 'String', matRad_playerInfoString(data, i, j, p, weight));
set(data.hTime, 'String', sprintf('%.1f / %.1f s (%s)', state.curTime, data.totalTime, data.timeKind));
set(data.hBeamValue, 'String', sprintf('%d / %d', i, data.numOfBeams));
set(data.hShapeValue, 'String', sprintf('%d / %d', j, max(nShapes, 1)));

% keep the sliders in sync with the animation, and re-range the shape
% slider whenever the beam (and with it its number of shapes) changes
if state.shownBeam ~= i
    matRad_playerSetSliderRange(data.hShapeSlider, nShapes);
end
matRad_playerSetSliderValue(data.hBeamSlider, i);
matRad_playerSetSliderValue(data.hShapeSlider, j);

state.shownBeam = i;
state.shownShape = j;
setappdata(hFig, 'matRadApertureState', state);
end

function str = matRad_playerInfoString(data, i, j, p, weight)

parts = {sprintf('w = %.4g', weight)};

if data.isVMAT
    if data.isDAOBeam(i)
        parts{end + 1} = 'DAO control point';
    else
        parts{end + 1} = 'interpolated control point';
    end
    mu = matRad_apertureFieldOrNaN(data.apertureInfo.beam(i).shape(j), 'MU');
    muRate = matRad_apertureFieldOrNaN(data.apertureInfo.beam(i).shape(j), 'MURate');
    if isfinite(mu)
        parts{end + 1} = sprintf('MU = %.4g', mu);
    end
    if isfinite(muRate)
        parts{end + 1} = sprintf('MU rate = %.4g MU/s', muRate);
    end
end

parts{end + 1} = sprintf('%s = %.3g s', data.timeLabel, data.duration(p));

str = strjoin(parts, '   |   ');
end

function matRad_playerSetSliderValue(hSlider, value)
value = min(max(value, get(hSlider, 'Min')), get(hSlider, 'Max'));
if get(hSlider, 'Value') ~= value
    set(hSlider, 'Value', value);
end
end

function matRad_playerSetSliderRange(hSlider, numSteps)
set(hSlider, 'Min', 1, 'Max', max(numSteps, 2), 'Value', 1);
if numSteps > 1
    step = 1 / (numSteps - 1);
    set(hSlider, 'SliderStep', [step, max(step, 0.1)], 'Enable', 'on');
else
    set(hSlider, 'Enable', 'off');
end
end

function matRad_playerSliderCb(hObj, ~, whichSlider)
% Scrubbing jumps playback to the start of the selected aperture, so that
% resuming continues from there.

hFig = ancestor(hObj, 'figure');
data = getappdata(hFig, 'matRadApertureData');
state = getappdata(hFig, 'matRadApertureState');

i = round(get(data.hBeamSlider, 'Value'));
j = round(get(data.hShapeSlider, 'Value'));
if strcmp(whichSlider, 'beam')
    % a new beam may carry fewer shapes than the one before
    j = min(j, max(matRad_apertureNumShapes(data.apertureInfo, i), 1));
end

p = find(data.panels(:, 1) == i & data.panels(:, 2) == j, 1);
if isempty(p)
    p = find(data.panels(:, 1) == i, 1);
end
if isempty(p)
    % beam without any shape - snap to the closest one that has one
    [~, p] = min(abs(data.panels(:, 1) - i));
end

state.curTime = data.tStart(p);
state.shownBeam = 0; % force the shape slider to be re-ranged
setappdata(hFig, 'matRadApertureState', state);

matRad_playerDraw(hFig);
end

function matRad_playerSpeedCb(hObj, ~)
hFig = ancestor(hObj, 'figure');
data = getappdata(hFig, 'matRadApertureData');
state = getappdata(hFig, 'matRadApertureState');
state.speed = data.speedList(get(hObj, 'Value'));
setappdata(hFig, 'matRadApertureState', state);
end

function matRad_playerToggleCb(hObj, ~)
hFig = ancestor(hObj, 'figure');
data = getappdata(hFig, 'matRadApertureData');
state = getappdata(hFig, 'matRadApertureState');
state.doLoop = logical(get(data.hLoop, 'Value'));
state.doInterp = logical(get(data.hInterp, 'Value'));
setappdata(hFig, 'matRadApertureState', state);
matRad_playerDraw(hFig);
end

function matRad_playerPlayCb(hObj, ~)
hFig = ancestor(hObj, 'figure');
state = getappdata(hFig, 'matRadApertureState');

if state.isPlaying
    % stops the loop below, which is still running in the callback that
    % started it and picks the flag up at its next drawnow
    state.isPlaying = false;
    setappdata(hFig, 'matRadApertureState', state);
else
    matRad_playerPlay(hFig);
end
end

function matRad_playerPlay(hFig)
% Advances the playback clock in real time until paused, until the plan
% ends (unless looping) or until the figure is closed.

state = getappdata(hFig, 'matRadApertureState');
if state.isPlaying
    return
end

data = getappdata(hFig, 'matRadApertureData');

% restart from the beginning when play is hit at the end of the plan
if state.curTime >= data.totalTime
    state.curTime = 0;
end
state.isPlaying = true;
setappdata(hFig, 'matRadApertureState', state);
set(data.hPlay, 'String', 'Pause');

tLast = tic;
while ishandle(hFig)

    state = getappdata(hFig, 'matRadApertureState');
    if ~state.isPlaying
        break
    end

    % integrate incrementally rather than from a fixed start, so that
    % changing the speed or scrubbing mid-flight is picked up right away
    dt = toc(tLast);
    tLast = tic;
    state.curTime = state.curTime + dt * state.speed;

    atEnd = state.curTime >= data.totalTime;
    if atEnd
        if state.doLoop
            state.curTime = mod(state.curTime, data.totalTime);
            atEnd = false;
        else
            state.curTime = data.totalTime - eps(data.totalTime);
            state.isPlaying = false;
        end
    end

    setappdata(hFig, 'matRadApertureState', state);
    matRad_playerDraw(hFig);
    drawnow;

    if atEnd
        break
    end

    pause(0.02); % ~50 fps at most, and yields to the UI
end

if ishandle(hFig)
    state = getappdata(hFig, 'matRadApertureState');
    state.isPlaying = false;
    setappdata(hFig, 'matRadApertureState', state);
    set(data.hPlay, 'String', 'Play');
end
end

%| pragma Justify (metric, "file_length",
%|                  "the aperture views and the player share the aperture " +
%|                  "geometry, the weight colouring and the panel " +
%|                  "bookkeeping; splitting the file up would only move " +
%|                  "those into helpers that nothing else uses");
