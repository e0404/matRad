function test_suite = test_visApertureInfo
% The output should always be test_suite, and the function name the same as
% your file name

% To collect all tests defined below, this is needed in newer Matlab
% versions. test_functions will collect function handles to below test
% functions
test_functions = localfunctions();

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

%% fixtures

function apertureInfo = helper_getStaticApertureInfo()
% A sequenced (static / step-and-shoot) plan, i.e. few beams carrying many
% segments each. Cached, since sequencing the test data is by far the most
% expensive part of this file and the visualization never modifies it.
persistent cachedInfo
if isempty(cachedInfo)
    p = load('photons_testData.mat');
    pln = p.pln;
    pln.propSeq.sequencer = 'siochi';
    resultGUI = matRad_sequencing(p.resultGUI, p.stf, pln);
    cachedInfo = resultGUI.sequencing.apertureInfo;
end
apertureInfo = cachedInfo;

function apertureInfo = helper_getVmatApertureInfo()
% A sequenced VMAT arc, i.e. many control points carrying one aperture each,
% part of them interpolated rather than optimized.
persistent cachedInfo
if isempty(cachedInfo)
    p = load('photons_testData.mat', 'ct', 'cst', 'pln');
    pln = p.pln;

    pln.propStf.generator                = 'PhotonVMAT';
    pln.propStf.gantryAngles             = [-180, 180];
    pln.propStf.couchAngles              = [0, 0];
    pln.propStf.maxGantryAngleSpacing    = 15;
    pln.propStf.maxDAOGantryAngleSpacing = 30;
    pln.propStf.maxFMOGantryAngleSpacing = 45;
    pln.propStf.isoCenter                = matRad_getIsoCenter(p.cst, p.ct, 0);

    pln.propSeq.continuousAperture = false;
    pln.propOpt.runVMAT = true;

    stf = matRad_generateStf(p.ct, p.cst, pln);

    sequencer = matRad_SequencingPhotonsSiochiLeaf(pln);
    sequencer.runVMAT = true;
    sequencer.numLevels = 5;
    sequencer.weightToMU = 100;

    w = ones(sum([stf.numOfRays]), 1);
    cachedInfo = sequencer.sequence(w, stf).apertureInfo;
end
apertureInfo = cachedInfo;

function apertureInfo = helper_getShortArc()
% A hand-built miniature arc. Playback runs in delivery time, so the tests
% that actually let the animation run need a plan that is over in well under
% a second rather than a realistic arc, which takes minutes to deliver.
numOfBeams = 5;
numOfLeafPairs = 4;

apertureInfo.bixelWidth = 5;
apertureInfo.runVMAT = true;

for i = 1:numOfBeams
    apertureInfo.beam(i).gantryAngle = mod(340 + (i - 1) * 30, 360); % wraps through 0
    apertureInfo.beam(i).MLCWindow = [-25 25];
    apertureInfo.beam(i).numOfActiveLeafPairs = numOfLeafPairs;
    apertureInfo.beam(i).isActiveLeafPair = ones(numOfLeafPairs, 1);
    apertureInfo.beam(i).leafPairPos = ((1:numOfLeafPairs)' - 2.5) * apertureInfo.bixelWidth;
    apertureInfo.beam(i).numOfShapes = 1;
    apertureInfo.beam(i).time = 0.05;
    apertureInfo.beam(i).gantryRot = 30 / 0.05;

    apertureInfo.beam(i).shape(1).leftLeafPos = -10 + i * ones(numOfLeafPairs, 1);
    apertureInfo.beam(i).shape(1).rightLeafPos = 10 + i * ones(numOfLeafPairs, 1);
    apertureInfo.beam(i).shape(1).weight = 0.2 * i;
    apertureInfo.beam(i).shape(1).MU = 20 * i;
    apertureInfo.beam(i).shape(1).MURate = 20 * i / 0.05;

    apertureInfo.arc.beam(i).isDAOBeam = mod(i, 2) == 1;
    apertureInfo.arc.beam(i).doseAngleBordersDiff = 30;
end

function figHandles = helper_openFigures()
figHandles = findobj(0, 'Type', 'figure');

function numNew = helper_closeNewFigures(figHandlesBefore)
% Number of figures opened since figHandlesBefore was taken; closes them, so
% that a test neither leaks figures nor closes any that it did not open.
newFigures = setdiff(findobj(0, 'Type', 'figure'), figHandlesBefore);
numNew = numel(newFigures);
close(newFigures);

function hFig = helper_theNewFigure(figHandlesBefore)
% The one figure opened since figHandlesBefore was taken. Deliberately not
% gcf: a test must not depend on which figure another test left current.
newFigures = setdiff(findobj(0, 'Type', 'figure'), figHandlesBefore);
assertEqual(numel(newFigures), 1, 'expected exactly one new figure');
hFig = newFigures;

function numTabs = helper_numTabs()
numTabs = numel(findobj(0, 'Type', 'uitab'));

function label = helper_yLabelOfFirstAxes(hFig)
allAxes = findobj(hFig, 'Type', 'axes');
assertFalse(isempty(allAxes), 'no axes were drawn');
label = get(get(allAxes(1), 'YLabel'), 'String');

%% every view has to open on a real plan

function test_staticViewsOpen
apertureInfo = helper_getStaticApertureInfo();

% 'auto' resolves to the per-beam layout for a static plan
before = helper_openFigures();
matRad_visApertureInfo(apertureInfo);
assertTrue(helper_closeNewFigures(before) > 0);

for view = {'perBeam', 'grid', 'interactive'}
    before = helper_openFigures();
    matRad_visApertureInfo(apertureInfo, view{1});
    assertTrue(helper_closeNewFigures(before) > 0, ...
               sprintf('view ''%s'' did not open a figure', view{1}));
end

function test_vmatViewsOpen
apertureInfo = helper_getVmatApertureInfo();
assertTrue(apertureInfo.runVMAT);

% 'auto' resolves to grid + trajectory + metrics for an arc
before = helper_openFigures();
matRad_visApertureInfo(apertureInfo);
assertTrue(helper_closeNewFigures(before) >= 3);

for view = {'perBeam', 'grid', 'trajectory', 'metrics', 'interactive'}
    before = helper_openFigures();
    matRad_visApertureInfo(apertureInfo, view{1});
    assertTrue(helper_closeNewFigures(before) > 0, ...
               sprintf('view ''%s'' did not open a figure', view{1}));
end

function test_bothLeafCoordinatesOpen
% The leaf coordinates change the vertical axis of every aperture plot: leaf
% pair index by default, physical position on request.
apertureInfo = helper_getStaticApertureInfo();

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'grid', 'leafCoordinate', 'leafNum');
assertEqual(helper_yLabelOfFirstAxes(helper_theNewFigure(before)), 'leaf pair #');
helper_closeNewFigures(before);

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'grid', 'leafCoordinate', 'physical');
assertEqual(helper_yLabelOfFirstAxes(helper_theNewFigure(before)), 'vert. pos. [mm]');
helper_closeNewFigures(before);

% option names and values are matched case-insensitively
before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'grid', 'LeafCoordinate', 'PHYSICAL');
assertEqual(helper_yLabelOfFirstAxes(helper_theNewFigure(before)), 'vert. pos. [mm]');
helper_closeNewFigures(before);

function test_severalViewsAtOnce
apertureInfo = helper_getVmatApertureInfo();

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, {'trajectory', 'metrics'});
assertEqual(helper_closeNewFigures(before), 2);

%% per-beam layout

function test_perBeamTabbedPerBeam
% One tab per beam in a single figure, or - where uitabgroup is missing, as
% in Octave - one figure per beam.
apertureInfo = helper_getStaticApertureInfo();
numOfBeams = numel(apertureInfo.beam);

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'perBeam');

if moxunit_util_platform_is_octave
    assertEqual(helper_closeNewFigures(before), numOfBeams);
else
    assertEqual(helper_numTabs(), numOfBeams);
    assertEqual(helper_closeNewFigures(before), 1);
end

function test_perBeamDrawsEveryShape
apertureInfo = helper_getStaticApertureInfo();

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'perBeam');

% every beam gets its shapes tiled, whether they live on a tab or a figure
if moxunit_util_platform_is_octave
    parents = setdiff(findobj(0, 'Type', 'figure'), before);
else
    parents = findobj(0, 'Type', 'uitab');
end
assertEqual(numel(parents), numel(apertureInfo.beam));

numOfAxes = sort(arrayfun(@(h) numel(findobj(h, 'Type', 'axes')), parents(:)));
expected = sort(arrayfun(@(b) numel(b.shape), apertureInfo.beam(:)));
assertEqual(numOfAxes, expected);

helper_closeNewFigures(before);

%% player

function test_interactiveBuildsTransport
% The player starts paused and offers the two sliders (control point and
% shape) plus the transport.
apertureInfo = helper_getVmatApertureInfo();

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'interactive');
hFig = helper_theNewFigure(before);

assertEqual(numel(findobj(hFig, 'Style', 'slider')), 2);
assertEqual(numel(findobj(hFig, 'Style', 'checkbox')), 2);
assertEqual(numel(findobj(hFig, 'Style', 'popupmenu')), 1);

hPlay = findobj(hFig, 'Style', 'pushbutton');
assertEqual(numel(hPlay), 1);
assertEqual(get(hPlay, 'String'), 'Play');

data = getappdata(hFig, 'matRadApertureData');
assertFalse(isempty(data));

% one entry per aperture, and a positive dwell for each of them
numOfShapes = sum(arrayfun(@(b) numel(b.shape), apertureInfo.beam));
assertEqual(size(data.panels, 1), numOfShapes);
assertEqual(numel(data.duration), numOfShapes);
assertTrue(all(isfinite(data.duration)) && all(data.duration > 0));
assertElementsAlmostEqual(data.totalTime, sum(data.duration));

% the control point slider spans the beams
hBeamSlider = findobj(hFig, 'Style', 'slider');
assertEqual(max([get(hBeamSlider(1), 'Max'), get(hBeamSlider(2), 'Max')]), ...
            numel(apertureInfo.beam));

assertEqual(helper_closeNewFigures(before), 1);

function test_playerUsesDeliveryTimeForVmat
% VMAT carries the machine timing, so the dwell of a control point is its
% delivery time rather than a weight-proportional stand-in.
before = helper_openFigures();
matRad_visApertureInfo(helper_getVmatApertureInfo(), 'interactive');
data = getappdata(helper_theNewFigure(before), 'matRadApertureData');
assertTrue(data.timeIsPhysical);
helper_closeNewFigures(before);

% a static plan has no timing, so the dwells fall back to the shape weights
before = helper_openFigures();
matRad_visApertureInfo(helper_getStaticApertureInfo(), 'interactive');
data = getappdata(helper_theNewFigure(before), 'matRadApertureData');
assertFalse(data.timeIsPhysical);
assertTrue(all(data.duration > 0));
helper_closeNewFigures(before);

function test_animatePlaysAndReturns
% 'animate' plays the plan once and hands control back, leaving the figure
% behind in the same state 'interactive' starts in.
apertureInfo = helper_getShortArc();

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'animate');
hFig = helper_theNewFigure(before);

data = getappdata(hFig, 'matRadApertureData');
state = getappdata(hFig, 'matRadApertureState');
assertFalse(state.isPlaying);
assertTrue(state.curTime > 0.9 * data.totalTime);
assertEqual(get(findobj(hFig, 'Style', 'pushbutton'), 'String'), 'Play');

assertEqual(helper_closeNewFigures(before), 1);

%% degenerate input

function test_zeroWeightPlanOpens
% A plan whose apertures all carry zero weight must not divide by zero,
% neither in the weight colouring nor in the playback timeline.
apertureInfo = helper_getShortArc();
for i = 1:numel(apertureInfo.beam)
    apertureInfo.beam(i).shape(1).weight = 0;
end
apertureInfo = rmfield(apertureInfo, 'arc');
apertureInfo.runVMAT = false;

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'interactive');
data = getappdata(helper_theNewFigure(before), 'matRadApertureData');
assertTrue(isfinite(data.totalTime) && data.totalTime > 0);
assertEqual(helper_closeNewFigures(before), 1);

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'grid');
assertEqual(helper_closeNewFigures(before), 1);

function test_vmatOnlyViewsSkippedForStaticPlan
% The trajectory and metrics views read apertureInfo.arc, which a static
% plan does not have - they warn and open nothing instead of erroring.
apertureInfo = helper_getStaticApertureInfo();

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'trajectory');
assertEqual(helper_closeNewFigures(before), 0);

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'metrics');
assertEqual(helper_closeNewFigures(before), 0);

%% argument handling

function test_invalidArgumentsError
apertureInfo = helper_getShortArc();

before = helper_openFigures();

% unknown view
assertExceptionThrown(@() matRad_visApertureInfo(apertureInfo, 'noSuchView'), 'matRad:Error');

% options have to come in name-value pairs
assertExceptionThrown(@() matRad_visApertureInfo(apertureInfo, 'grid', 'leafCoordinate'), ...
                      'matRad:Error');

% unknown option name
assertExceptionThrown(@() matRad_visApertureInfo(apertureInfo, 'grid', 'noSuchOption', 1), ...
                      'matRad:Error');

% unknown option value
assertExceptionThrown(@() matRad_visApertureInfo(apertureInfo, 'grid', 'leafCoordinate', 'nonsense'), ...
                      'matRad:Error');

helper_closeNewFigures(before);

function test_legacyCallSupported
% Deprecated signature matRad_visApertureInfo(apertureInfo,leafCoord,view):
% the leaf coordinates used to be the second argument, which is now the view.
apertureInfo = helper_getStaticApertureInfo();

before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'physical', 'grid');
assertEqual(helper_yLabelOfFirstAxes(helper_theNewFigure(before)), 'vert. pos. [mm]');
assertEqual(helper_closeNewFigures(before), 1);

% without the third argument the view falls back to 'auto', which for a
% static plan is the per-beam layout - a single tabbed figure, or one figure
% per beam where tabs are unavailable
before = helper_openFigures();
matRad_visApertureInfo(apertureInfo, 'physical');
newFigures = setdiff(findobj(0, 'Type', 'figure'), before);
assertEqual(helper_yLabelOfFirstAxes(newFigures(1)), 'vert. pos. [mm]');
assertTrue(helper_closeNewFigures(before) > 0);
