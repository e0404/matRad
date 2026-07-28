function matRad_visAperturePlayer(apertureInfo, leafCoord, wMax, isVMAT, numOfBeams, autoPlay)
% matRad_visAperturePlayer Single-aperture player of matRad_visApertureInfo.
%
%   One axes showing one aperture at a time, driven either by the transport
%   (play/pause, speed, loop) or by the two sliders (control point / beam,
%   and shape within that beam). The 'animate' and 'interactive' views are
%   the same figure, they only differ in whether it starts playing - pausing
%   an animation leaves the sliders free to scrub.
%
% call:
%   matRad_visAperturePlayer(apertureInfo, leafCoord, wMax, isVMAT, numOfBeams, autoPlay)
%
% input:
%   apertureInfo: aperture weight and shape info struct
%   leafCoord:    'leafNum' or 'physical', see matRad_visApertureInfo
%   wMax:         weight the colour scale saturates at
%   isVMAT:       whether the plan is delivered as an arc
%   numOfBeams:   number of beams of the plan
%   autoPlay:     start playing ('animate') rather than paused ('interactive')
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

    % the leaf rows of a beam do not depend on the shape drawn on it, so any
    % shape gives the vertical extent
    [~, ~, yRow] = matRad_apertureQuads(apertureInfo, i, ...
                                        apertureInfo.beam(i).shape(1).leftLeafPos, ...
                                        apertureInfo.beam(i).shape(1).rightLeafPos, ...
                                        leafCoord, bixelWidth);
    yLo = min(yLo, min(yRow(:)));
    yHi = max(yHi, max(yRow(:)));
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
