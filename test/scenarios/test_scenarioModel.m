function test_suite = test_scenarioModel

test_functions=localfunctions();

initTestSuite;

%Add automated instance tests for avaiable models
models = {matRad_NominalScenario(),matRad_WorstCaseScenarios(),matRad_ImportanceScenarios(),matRad_RandomScenarios()};

instanceTests = cellfun(@func2str,test_functions,'UniformOutput',false);
funIx = ~cellfun(@isempty,strfind(instanceTests,'instanceTest_'));
instanceTests = instanceTests(funIx);
instanceTestHandles = cellfun(@str2func,instanceTests,'UniformOutput',false);

for m = 1:numel(models)
    model = models{m};
    for i = 1:numel(instanceTests)
        funHandle = @() instanceTestHandles{i}(model);
        test_case=MOxUnitFunctionHandleTestCase([instanceTests{i} '_' class(model)],mfilename,funHandle);
        test_suite=addTest(test_suite, test_case); 
    end           
end

function test_scenarioAbstract
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() matRad_ScenarioModel(),'');
    else
        assertExceptionThrown(@() matRad_ScenarioModel(),'MATLAB:class:abstract');
    end

function test_scenarioAbstractAvailableTypes()
    availableTypes = matRad_ScenarioModel.AvailableScenCreationTYPE();
    assertTrue(iscell(availableTypes))
    assertTrue(all(cellfun(@ischar,availableTypes)));

    for i = 1:numel(availableTypes)
        model = matRad_multScen([],availableTypes{i});
        assertTrue(isa(model,'matRad_ScenarioModel'));
        assertEqual(model.name,availableTypes{i});
    end


function instanceTest_listAllScenarios(model)
    model.listAllScenarios();
    assertTrue(true); %Will be reached if above call does not fail

function instanceTest_relRangeUncertainty(model)
    newValue = 0.01;
    model.rangeRelSD = newValue;
    getValue = model.rangeRelSD;
    assertEqual(newValue,getValue);   

function instanceTest_absRangeUncertainty(model)
    newValue = 2;
    model.rangeAbsSD = newValue;
    getValue = model.rangeAbsSD;
    assertEqual(newValue,getValue); 

function instanceTest_shiftUncertainty(model)
    newValue = [1 1 1];
    model.shiftSD = newValue;
    getValue = model.shiftSD;
    assertEqual(newValue,getValue); 

function instanceTest_wcSigma(model)
    newValue = 5;
    model.wcSigma = newValue;
    getValue = model.wcSigma;
    assertEqual(newValue,getValue); 

function instanceTest_phaseProbability(model)
    newValue = rand(model.numOfCtScen,1);
    newValue = newValue ./ sum(newValue);
    model.phaseProbability = newValue;
    getValue = model.phaseProbability;
    assertEqual(newValue,getValue); 

function instanceTest_TYPE(model)
    %assertWarning(@() model.TYPE,'matRad:Deprecated');
    assertEqual(model.TYPE,model.name);

function instanceTest_wcFactor(model)
    %assertWarning(@() model.wcFactor,'matRad:Deprecated');
    assertEqual(model.TYPE,model.name);

