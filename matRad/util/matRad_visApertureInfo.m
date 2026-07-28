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
