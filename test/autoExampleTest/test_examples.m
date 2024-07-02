function test_suite = test_examples

matRad_cfg = MatRad_Config.instance();

% supressing the inherent Ocatave warnings for division by zero
if matRad_cfg.isOctave
    warning('off','Octave:divide-by-zero');
end

% Define Scripts relative to root folder

exampleScripts = {'examples/matRad_example1_phantom.m',...
    'examples/matRad_example2_photons.m',...
    'examples/matRad_example3_photonsDAO.m',...
    'examples/matRad_example4_photonsMC.m',...
    'examples/matRad_example5_protons.m',...
    'examples/matRad_example6_protonsNoise.m',...
    'examples/matRad_example7_carbon.m',... 
    'examples/matRad_example8_protonsRobust.m',...
    'examples/matRad_example9_4DDoseCalcMinimal.m',... 
    'examples/matRad_example10_4DphotonRobust.m',...
    'examples/matRad_example11_helium.m',...
    'examples/matRad_example12_simpleParticleMonteCarlo.m',...
    'matRad.m',...
    };

%exampleScripts = {'matRad.m'}; %Uncomment to fast-test the example testing workflow

exampleScripts = cellfun(@(script) fullfile(matRad_cfg.matRadRoot,script),exampleScripts,'UniformOutput',false);


testing_prefix = 'tmptest_';

% Some parameters to reduce computational overhead during testing
unitTestBixelWidth = 20;
unitTestSpotSpacing = matRad_cfg.defaults.propStf.longitudinalSpotSpacing;
unitTestResolution = matRad_cfg.defaults.propDoseCalc.doseGrid.resolution;

%% Copy and manipulate all scripts
[folders,names,exts] = cellfun(@fileparts,exampleScripts,'UniformOutput',false);

%Create temporary example test folder
tmpExampleTestFolder = tempdir();
tmpExampleTestFolder = fullfile(tmpExampleTestFolder,'exampleTest');
if ~exist(tmpExampleTestFolder,'dir')
    mkdir(tmpExampleTestFolder);
end
addpath(tmpExampleTestFolder);
newFolders = cell(size(folders));
[newFolders{:}] = deal(tmpExampleTestFolder);

%Copy scripts
testScriptNames = strcat(testing_prefix,names);  
testScriptFiles = strcat(testScriptNames,exts);
testScripts = cellfun(@fullfile,newFolders,testScriptFiles,'UniformOutput',false);

status = cellfun(@copyfile,exampleScripts,testScripts);

matRad_unitTestTextManipulation(testScriptFiles,'pln.propStf.bixelWidth',['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth) ';'],tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles,'pln.propStf.longitudinalSpotSpacing',['pln.propStf.longitudinalSpotSpacing = ' num2str(unitTestSpotSpacing) ';'],tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles,'pln.propDoseCalc.doseGrid.resolution.x',['pln.propDoseCalc.doseGrid.resolution.x = ' num2str(unitTestResolution.x) ';'],tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles,'pln.propDoseCalc.doseGrid.resolution.y',['pln.propDoseCalc.doseGrid.resolution.y = ' num2str(unitTestResolution.y) ';'],tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles,'pln.propDoseCalc.doseGrid.resolution.z',['pln.propDoseCalc.doseGrid.resolution.z = ' num2str(unitTestResolution.z) ';'],tmpExampleTestFolder);
matRad_unitTestTextManipulation(testScriptFiles,'display(','%%%%%%%%%%%%%%% REMOVED DISPLAY FOR TESTING %%%%%%%%%%%%%%',tmpExampleTestFolder);

%initTestSuite;
%We need to manually set up the test_suite to bypass the automatic function
%assignment
test_suite=MOxUnitTestSuite();

for testIx = 1:length(testScriptNames)
    currScriptName = testScriptNames{testIx};
    testfun = @() runSingleExampleTest(currScriptName,testScripts{testIx}); %Test is evaluated in the base workspace and clears new variables after that

    test_case=MOxUnitFunctionHandleTestCase(...
        names{testIx},...
        mfilename, testfun);
    test_suite=addTest(test_suite, test_case);
    %test_functions{testIx,1} = testfun;
end

try
    rmdir(exampleTestFolder,'s');
catch
    warning('Could not delete temporary example test folder');
end
    
%initTestSuite;
%We need to manually set up the test_suite

function runSingleExampleTest(exampleName,path)
%We use a kind of fishy way to run the test in the base workspace
%First we record the variables we have in the base workspace
baseVars = evalin('base','who');

%Example is evaluated in the base workspace
%addpath(fileparts(path));
try
    evalin('base',exampleName);

    %Clean up of the base workspace by cleaning all new variables
    afterTestVars = evalin('base','who');
    newVars = setdiff(afterTestVars,baseVars);
    evalin('base',['clear ' strjoin(newVars)]);
    close all;
catch ME
    %Also clean up the base workspace by cleaning all new variables
    afterTestVars = evalin('base','who');
    newVars = setdiff(afterTestVars,baseVars);
    evalin('base',['clear ' strjoin(newVars)]);
    close all;

    %Now rethrow
    rethrow(ME);
end