function test_suite = test_examples

matRad_cfg = MatRad_Config.instance();
matRad_cfg.setDefaultPropertiesForTesting();

% supressing the inherent Ocatave warnings for division by zero
if matRad_cfg.isOctave
    warning('off','Octave:divide-by-zero');
end

% Define Scripts
exampleScripts = {'../examples/matRad_example1_phantom.m',...
    '../examples/matRad_example2_photons.m',...
    '../examples/matRad_example3_photonsDAO.m',...
    '../examples/matRad_example4_photonsMC.m',...
    '../examples/matRad_example5_protons.m',...
    '../examples/matRad_example6_protonsNoise.m',...
    '../examples/matRad_example7_carbon.m',... 
    '../examples/matRad_example8_protonsRobust.m',...
    '../examples/matRad_example9_4DDoseCalcMinimal.m',... 
    '../examples/matRad_example10_4DphotonRobust.m',...
    '../examples/matRad_example11_helium.m',...
    '../matRad.m',...
    };

testing_prefix = 'test_';

% Some parameters to reduce computational overhead during testing
unitTestBixelWidth = 20;
unitTestSpotSpacing = matRad_cfg.propStf.defaultLongitudinalSpotSpacing;
unitTestResolution = matRad_cfg.propDoseCalc.defaultResolution;

% Copy and manipulate all scripts
[folders,names,exts] = cellfun(@fileparts,exampleScripts,'UniformOutput',false);
testScriptNames = strcat(testing_prefix,names);    
testScripts = cellfun(@fullfile,folders,strcat(testScriptNames,exts),'UniformOutput',false);
status = cellfun(@copyfile,exampleScripts,testScripts);

matRad_unitTestTextManipulation(testScripts,'pln.propStf.bixelWidth',['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth)]);
matRad_unitTestTextManipulation(testScripts,'pln.propStf.longitudinalSpotSpacing',['pln.propStf.longitudinalSpotSpacing = ' num2str(unitTestBixelWidth)]);
matRad_unitTestTextManipulation(testScripts,'pln.propDoseCalc.resolution.x',['pln.propDoseCalc.resolution.x = ' num2str(unitTestResolution.x)]);
matRad_unitTestTextManipulation(testScripts,'pln.propDoseCalc.resolution.y',['pln.propDoseCalc.resolution.y = ' num2str(unitTestResolution.y)]);
matRad_unitTestTextManipulation(testScripts,'pln.propDoseCalc.resolution.z',['pln.propDoseCalc.resolution.z = ' num2str(unitTestResolution.z)]);
matRad_unitTestTextManipulation(testScripts,'display(','%%%%%%%%%%%%%%% REMOVED DISPLAY FOR TESTING %%%%%%%%%%%%%%');

%initTestSuite;
%We need to manually set up the test_suite to bypass the automatic function
%assignment
test_suite=MOxUnitTestSuite();

for testIx = 1:length(testScriptNames)
    currScriptName = testScriptNames{testIx};
    testfun = @() evalin('caller',currScriptName); %this is ridiculous style but it works. We trust that we don't have a variable in the workspace that overwrites the callers (MOxUnit) variable
    test_case=MOxUnitFunctionHandleTestCase(...
        names{testIx},...
        mfilename, testfun);
    test_suite=addTest(test_suite, test_case);
    %test_functions{testIx,1} = testfun;
end
    
%initTestSuite;
%We need to manually set up the test_suite





