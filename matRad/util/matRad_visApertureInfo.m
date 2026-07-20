function matRad_visApertureInfo(apertureInfo, mode, viewType)
% matRad function to visualize aperture shapes stored as struct
%
% call:
%   matRad_visApertureInfo(apertureInfo)
%   matRad_visApertureInfo(apertureInfo,mode)
%   matRad_visApertureInfo(apertureInfo,mode,viewType)
%
% input:
%   apertureInfo: aperture weight and shape info struct
%   mode:         switch to display leaf numbers ('leafNum') or physical
%                 coordinates of the leaves ('physical'). Default: 'leafNum'
%   viewType:         layout to use. Default: 'auto'
%                 'auto'       - 'perBeam' for static (IMRT/DAO) plans,
%                                'grid' + 'trajectory' + 'metrics' for VMAT
%                 'perBeam'    - one figure per beam, one subplot per shape.
%                                Suited to static delivery, where a beam
%                                carries many segments
%                 'grid'       - one (paginated) figure, one subplot per
%                                control point. Suited to VMAT, where every
%                                control point carries a single aperture
%                 'trajectory' - leaf position vs gantry angle (VMAT only)
%                 'metrics'    - MU rate, gantry rotation speed and leaf
%                                speed vs gantry angle, against the machine
%                                delivery constraints (VMAT only)
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
    mode = 'leafNum'; % options: 'physical','leafNum'
end
if nargin < 3 || isempty(viewType)
    viewType = 'auto';
end

isVMAT = isfield(apertureInfo, 'runVMAT') && apertureInfo.runVMAT;

% resolve 'auto': the per-beam layout tiles the shapes of a single beam,
% which is the right thing for static delivery (few beams, many segments
% each) but degenerates for VMAT (many control points, one aperture each -
% one figure per control point). VMAT therefore gets the control-point grid
% plus the two views that actually show the arc: leaf trajectories and the
% delivery metrics against the machine constraints.
if strcmp(viewType, 'auto')
    if isVMAT
        viewList = {'grid', 'trajectory', 'metrics'};
    else
        viewList = {'perBeam'};
    end
else
    viewList = {viewType};
end

% global parameters
numOfBeams = size(apertureInfo.beam, 2);
wMax = matRad_getMaxWeight(apertureInfo, isVMAT);

for v = 1:numel(viewList)
    switch viewList{v}
        case 'perBeam'
            matRad_visPerBeam(apertureInfo, mode, wMax, isVMAT, numOfBeams);
        case 'grid'
            matRad_visGrid(apertureInfo, mode, wMax, numOfBeams);
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

function color = matRad_weightColormap()
% custom colormap: light -> saturated red with increasing weight
color = [0.2:0.01:0.8; 0.2:0.01:0.8; 0.2:0.01:0.8]';
color = flipud(color);
color(:, 3) = 0;
color(:, 2) = 0;
end

function n = matRad_numShapes(apertureInfo, i)
% Deliberately count the shape structs rather than reading
% beam(i).numOfShapes: in VMAT, interpolated (non-DAO) beams carry a
% computed shape while their numOfShapes (the number of *optimized* shapes)
% is 0 - and those interpolated apertures should be shown too.
n = numel(apertureInfo.beam(i).shape);
end

function matRad_drawAperture(apertureInfo, i, j, mode, bixelWidth)
% Draws the closed leaf pairs of shape j of beam i into the current axes.

minX = apertureInfo.beam(i).MLCWindow(1);
maxX = apertureInfo.beam(i).MLCWindow(2);

if strcmp(mode, 'leafNum')
    % the leaf indices have to be flipped in order to fit to the order of
    % the leaf positions (1st row of leafPos is lowest row in physical
    % coordinates)
    activeLeafInd = flipud(find(apertureInfo.beam(i).isActiveLeafPair));
end

for k = 1:apertureInfo.beam(i).numOfActiveLeafPairs

    leftPos = apertureInfo.beam(i).shape(j).leftLeafPos(k);
    rightPos = apertureInfo.beam(i).shape(j).rightLeafPos(k);

    if strcmp(mode, 'physical')
        rowCentre = apertureInfo.beam(i).leafPairPos(k);
        halfRow = bixelWidth / 2;
    else
        rowCentre = activeLeafInd(k);
        halfRow = 1 / 2;
    end

    rowLo = rowCentre - halfRow;
    rowHi = rowCentre + halfRow;

    % left leaf, then right leaf
    fill([minX leftPos leftPos minX], [rowLo rowLo rowHi rowHi], [0.5 0.5 0.5]);
    fill([rightPos maxX maxX rightPos], [rowLo rowLo rowHi rowHi], [0.5 0.5 0.5]);
end

axis tight;
xlabel('horiz. pos. [mm]');
if strcmp(mode, 'physical')
    ylabel('vert. pos. [mm]');
else
    ylabel('leaf pair #');
end
end

function matRad_visPerBeam(apertureInfo, mode, wMaxGlobal, isVMAT, numOfBeams)
% Legacy layout: one figure per beam, one subplot per shape of that beam.

bixelWidth = apertureInfo.bixelWidth;
color = matRad_weightColormap();

for i = 1:numOfBeams

    nShapes = matRad_numShapes(apertureInfo, i);

    % open new figure for every beam
    figure('units', 'inches');

    if isVMAT
        wMax = wMaxGlobal;
    else
        % if not VMAT, let wMax be the max weight of a particular angle
        wMax = wMaxGlobal;
        if nShapes > 0
            wMax = max([apertureInfo.beam(i).shape(:).weight]);
        end
    end

    subplotColumns = max(ceil(nShapes / 2), 1);
    subplotLines = max(ceil(nShapes / subplotColumns), 1);

    % adjust figure position
    set(gcf, 'pos', [0 0 1.8 * subplotColumns 3 * subplotLines]);

    % loop over all shapes of the beam
    for j = 1:nShapes

        subplot(subplotLines, subplotColumns, j);

        title(['Beam: ' num2str(i) ' Shape: ' num2str(j) ' w=' ...
               num2str(apertureInfo.beam(i).shape(j).weight, 2)], ...
              'Fontsize', 8);

        set(gca, 'Color', matRad_shapeColor(color, apertureInfo.beam(i).shape(j).weight, wMax));
        hold on;

        matRad_drawAperture(apertureInfo, i, j, mode, bixelWidth);
    end
end
end

function matRad_visGrid(apertureInfo, mode, wMax, numOfBeams)
% One subplot per control point, paginated over as many figures as needed.
% This is the VMAT counterpart of matRad_visPerBeam: there the shapes of one beam
% are tiled, here the (single) aperture of every beam is.

bixelWidth = apertureInfo.bixelWidth;
color = matRad_weightColormap();
maxPerFigure = 30; % a 6x5 grid is still legible, and fits a typical single arc

% flatten all (beam, shape) pairs into one list of panels
panels = zeros(0, 2);
for i = 1:numOfBeams
    for j = 1:matRad_numShapes(apertureInfo, i)
        panels(end + 1, :) = [i j]; %#ok<AGROW>
    end
end

numPanels = size(panels, 1);
if numPanels == 0
    return
end

% Spread the panels evenly over the pages instead of filling each to the
% brim, so that e.g. 31 apertures give two balanced pages rather than a
% full one followed by a page holding a single aperture.
numPages = ceil(numPanels / maxPerFigure);
perPage = ceil(numPanels / numPages);
isDAOBeam = matRad_getDAOFlags(apertureInfo, numOfBeams);

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

        subplot(subplotLines, subplotColumns, p - first + 1);

        % Label by gantry angle rather than beam index - along an arc the
        % angle is what identifies a control point. Interpolated (non-DAO)
        % control points are marked, since their apertures are derived
        % rather than optimized.
        angleLabel = sprintf('G=%.0f%s', apertureInfo.beam(i).gantryAngle, char(176));
        if ~isDAOBeam(i)
            angleLabel = [angleLabel ' (interp)'];
        end
        title(sprintf('%s w=%s', angleLabel, ...
                      num2str(apertureInfo.beam(i).shape(j).weight, 2)), 'Fontsize', 7);

        set(gca, 'Color', matRad_shapeColor(color, apertureInfo.beam(i).shape(j).weight, wMax));
        hold on;

        matRad_drawAperture(apertureInfo, i, j, mode, bixelWidth);
        set(gca, 'Fontsize', 6);
    end
end
end

function c = matRad_shapeColor(color, weight, wMax)
colorInd = max(ceil((weight / wMax) * (size(color, 1) - 1) + eps), 1);
colorInd = min(colorInd, size(color, 1));
c = color(colorInd, :);
end

function isDAOBeam = matRad_getDAOFlags(apertureInfo, numOfBeams)
isDAOBeam = true(1, numOfBeams);
if isfield(apertureInfo, 'arc') && isfield(apertureInfo.arc, 'beam')
    isDAOBeam = logical([apertureInfo.arc.beam.isDAOBeam]);
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
    if matRad_numShapes(apertureInfo, i) < 1
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
isDAOBeam = matRad_getDAOFlags(apertureInfo, numOfBeams);
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

isDAOBeam = matRad_getDAOFlags(apertureInfo, numOfBeams);

hasShape = false(1, numOfBeams);
for i = 1:numOfBeams
    hasShape(i) = matRad_numShapes(apertureInfo, i) >= 1;
end
beamIx = find(isDAOBeam & hasShape);

angles = zeros(1, numel(beamIx));
muRate = zeros(1, numel(beamIx));
gantryRot = zeros(1, numel(beamIx));
leafSpeed = zeros(1, numel(beamIx));

for b = 1:numel(beamIx)
    i = beamIx(b);
    angles(b) = apertureInfo.beam(i).gantryAngle;
    muRate(b) = matRad_scalarOrNaN(apertureInfo.beam(i).shape(1), 'MURate');
    gantryRot(b) = matRad_scalarOrNaN(apertureInfo.beam(i), 'gantryRot');
    leafSpeed(b) = matRad_scalarOrNaN(apertureInfo.beam(i), 'maxLeafSpeed');
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

function val = matRad_scalarOrNaN(s, fieldName)
val = NaN;
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName)(1);
end
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
