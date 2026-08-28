function test_suite = test_examples

matRad_cfg = MatRad_Config.instance();

% suppressing the inherent Ocatave warnings for division by zero
if matRad_cfg.isOctave
    warning('off', 'Octave:divide-by-zero');
end

% Define Scripts relative to root folder

exampleScripts = {'examples/matRad_example1_phantom.m', ...
                  'examples/matRad_example2_photons.m', ...
                  'examples/matRad_example3_photonsDAO.m', ...
                  'examples/matRad_example4_photonsMC.m', ...
                  'examples/matRad_example5_protons.m', ...
                  'examples/matRad_example6_protonsNoise.m', ...
                  'examples/matRad_example7_carbon.m', ...
                  'examples/matRad_example8_protonsRobust.m', ...
                  'examples/matRad_example9_4DDoseCalcMinimal.m', ...
                  'examples/matRad_example10_4DphotonRobust.m', ...
                  'examples/matRad_example11_helium.m', ...
                  'examples/matRad_example12_simpleParticleMonteCarlo.m', ...
                  'examples/matRad_example15_brachy.m', ...
                  'examples/matRad_example17_biologicalModels.m', ...
                  'examples/matRad_example19_CT_sCT_DVH_difference_photons.m', ...
                  'examples/matRad_example20_VHEE.m', ...
                  'examples/matRad_example22_photonsVMAT.m', ...
                  'examples/matRad_example23_neutrons.m', ...
                  'matRad.m' ...
                 };

% exampleScripts = {'matRad.m'}; %Uncomment to fast-test the example testing workflow

exampleScripts = cellfun(@(script) fullfile(matRad_cfg.matRadRoot, script), exampleScripts, 'UniformOutput', false);

testing_prefix = 'tmptest_';

% Some parameters to reduce computational overhead during testing
unitTestBixelWidth = 20;
unitTestSpotSpacing = matRad_cfg.defaults.propStf.longitudinalSpotSpacing;
unitTestResolution = matRad_cfg.defaults.propDoseCalc.doseGrid.resolution;

% Arc (VMAT) examples are by far the most expensive ones, since the number of
% gantry angles multiplies the cost of dose calculation, fluence optimization
% and direct aperture optimization alike. Collapse the arc to a handful of
% angles: one aperture per fluence-optimized angle and no interpolated control
% points in between. Every code path is still taken, only with a few beams.
unitTestArcSpacing = 120;               % [deg] -> 3 gantry angles on a full arc
unitTestAperturesPerFMOBeam = 1;        % one DAO control point per FMO angle
unitTestRefinedArcSpacing = 60;         % [deg] for the fine-arc recalculation
unitTestApertureVisMode = 'perBeam';    % static view; 'animate' plays the arc back in real time

%% Copy and manipulate all scripts
[folders, names, exts] = cellfun(@fileparts, exampleScripts, 'UniformOutput', false);

% Create temporary example test folder
tmpExampleTestFolder = helper_temporaryFolder('exampleTest', true);
addpath(tmpExampleTestFolder);
newFolders = cell(size(folders));
[newFolders{:}] = deal(tmpExampleTestFolder);

% Copy scripts
testScriptNames = strcat(testing_prefix, names);
testScriptFiles = strcat(testScriptNames, exts);
testScripts = cellfun(@fullfile, newFolders, testScriptFiles, 'UniformOutput', false);

status = cellfun(@copyfile, exampleScripts, testScripts);

matRad_unitTestTextManipulation(testScriptFiles, 'pln.propStf.bixelWidth', ...
                                ['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propStf.longitudinalSpotSpacing', ...
                                ['pln.propStf.longitudinalSpotSpacing = ' num2str(unitTestSpotSpacing) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propDoseCalc.doseGrid.resolution.x', ...
                                ['pln.propDoseCalc.doseGrid.resolution.x = ' num2str(unitTestResolution.x) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propDoseCalc.doseGrid.resolution.y', ...
                                ['pln.propDoseCalc.doseGrid.resolution.y = ' num2str(unitTestResolution.y) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propDoseCalc.doseGrid.resolution.z', ...
                                ['pln.propDoseCalc.doseGrid.resolution.z = ' num2str(unitTestResolution.z) ';'], tmpExampleTestFolder);

% Arc geometry. The three spacings must satisfy dose <= DAO <= FMO, and are all
% set to the same value here so that every dose-calculation angle is also a DAO
% control point and an FMO angle - no interpolated beams, one aperture each.
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propStf.maxGantryAngleSpacing', ...
                                ['pln.propStf.maxGantryAngleSpacing = ' num2str(unitTestArcSpacing) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propStf.maxDAOGantryAngleSpacing', ...
                                ['pln.propStf.maxDAOGantryAngleSpacing = ' num2str(unitTestArcSpacing) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propStf.maxFMOGantryAngleSpacing', ...
                                ['pln.propStf.maxFMOGantryAngleSpacing = ' num2str(unitTestArcSpacing) ';'], tmpExampleTestFolder);
% without this the stf generator subdivides the DAO angles again to reach the
% requested number of apertures per FMO beam, undoing the coarsening above
matRad_unitTestTextManipulation(testScriptFiles, 'pln.propStf.minAperturesPerFMOBeam', ...
                                ['pln.propStf.minAperturesPerFMOBeam = ' num2str(unitTestAperturesPerFMOBeam) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'refinedSpacing =', ...
                                ['refinedSpacing = ' num2str(unitTestRefinedArcSpacing) ';'], tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles, 'visMode = ''', ...
                                ['visMode = ''' unitTestApertureVisMode ''';'], tmpExampleTestFolder);

matRad_unitTestTextManipulation(testScriptFiles, 'disp(', ...
                                '%%%%%%%%%%%%%%% REMOVED DISPLAY FOR TESTING %%%%%%%%%%%%%%', tmpExampleTestFolder);

% initTestSuite;
% We need to manually set up the test_suite to bypass the automatic function
% assignment
test_suite = MOxUnitTestSuite();

for testIx = 1:length(testScriptNames)
    currScriptName = testScriptNames{testIx};
    % Test is evaluated in the base workspace and clears new variables after that
    testfun = @() helper_runSingleExampleTest(currScriptName, testScripts{testIx});

    test_case = MOxUnitFunctionHandleTestCase( ...
                                              names{testIx}, ...
                                              mfilename, testfun);
    test_suite = addTest(test_suite, test_case);
    % test_functions{testIx,1} = testfun;
end

% initTestSuite;
% We need to manually set up the test_suite

function helper_runSingleExampleTest(exampleName, path)
% We use a kind of fishy way to run the test in the base workspace
% First we record the variables we have in the base workspace
baseVars = evalin('base', 'who');

% Example is evaluated in the base workspace
% addpath(fileparts(path));
try
    evalin('base', exampleName);

    % Clean up of the base workspace by cleaning all new variables
    afterTestVars = evalin('base', 'who');
    newVars = setdiff(afterTestVars, baseVars);
    evalin('base', ['clear ' strjoin(newVars)]);
    close all;
catch ME
    % Also clean up the base workspace by cleaning all new variables
    afterTestVars = evalin('base', 'who');
    newVars = setdiff(afterTestVars, baseVars);
    evalin('base', ['clear ' strjoin(newVars)]);
    close all;

    % Now rethrow
    rethrow(ME);
end
