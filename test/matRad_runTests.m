%% This file runs the complete matRad test suite.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set path
run(['..' filesep 'matRad_rc']);

%% Prepare settings for testing
matRad_cfg = MatRad_Config.instance();
matRad_cfg.setDefaultPropertiesForTesting();

% supressing the inherent Ocatave warnings for division by zero
if strcmp(matRad_getEnvironment,'OCTAVE')
    warning('off','Octave:divide-by-zero');
end

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

testing_suffix = '_test';
unitTestBixelWidth = 20;
unitTestSpotSpacing = matRad_cfg.propStf.defaultLongitudinalSpotSpacing;
unitTestResolution = matRad_cfg.propDoseCalc.defaultResolution;

%Copy and manipulate all scripts
[folders,names,exts] = cellfun(@fileparts,exampleScripts,'UniformOutput',false);
testScriptNames = strcat(names,testing_suffix);    
testScripts = cellfun(@fullfile,folders,strcat(testScriptNames,exts),'UniformOutput',false);
status = cellfun(@copyfile,exampleScripts,testScripts);

matRad_unitTestTextManipulation(testScripts,'pln.propStf.bixelWidth',['pln.propStf.bixelWidth = ' num2str(unitTestBixelWidth)]);
matRad_unitTestTextManipulation(testScripts,'pln.propStf.longitudinalSpotSpacing',['pln.propStf.longitudinalSpotSpacing = ' num2str(unitTestBixelWidth)]);
matRad_unitTestTextManipulation(testScripts,'pln.propDoseCalc.resolution.x',['pln.propDoseCalc.resolution.x = ' num2str(unitTestResolution.x)]);
matRad_unitTestTextManipulation(testScripts,'pln.propDoseCalc.resolution.y',['pln.propDoseCalc.resolution.y = ' num2str(unitTestResolution.y)]);
matRad_unitTestTextManipulation(testScripts,'pln.propDoseCalc.resolution.z',['pln.propDoseCalc.resolution.z = ' num2str(unitTestResolution.z)]);
matRad_unitTestTextManipulation(testScripts,'display(','%%%%%%%%%%%%%%% REMOVED DISPLAY FOR TESTING %%%%%%%%%%%%%%');

errors = {};
%Running tests
for testIx = 1:length(testScriptNames)
    fprintf('Running Integration Test for ''%s''\n',names{testIx});
    try
        run(testScripts{testIx});
        clear ct cst pln stf dij resultGUI; %Make sure the workspace is somewhat clean
        delete(testScripts{testIx}); %Delete after successful run
    catch ME
        [~,scriptName] = fileparts(testScripts{testIx});
        if matRad_cfg.isMatlab
            message = ME.getReport();
        else
            message = ME.message;
        end
        errMsg = sprintf('Experiencd an error during testing of %s. Error-Message:\n %s',scriptName,message);
        warning(errMsg);
        errors{end+1} = errMsg;
    end
end

%Check if at least one script failed and report error
if ~isempty(errors)
    error(strjoin(errors,'\n\n============================\n\n'));
end
    