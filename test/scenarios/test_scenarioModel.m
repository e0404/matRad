function test_suite = test_scenarioModel

test_functions=localfunctions();

initTestSuite;

%Add automated instance tests for avaiable models
ct.numOfCtScen = 5;
models = {matRad_NominalScenario(),matRad_WorstCaseScenarios(),matRad_ImportanceScenarios(),matRad_RandomScenarios(),...
    matRad_NominalScenario(ct),matRad_WorstCaseScenarios(ct),matRad_ImportanceScenarios(ct),matRad_RandomScenarios(ct)};

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

function assignmentTestHelper(model,property,value)
    model.(property) = value;

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
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeRelSD','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeRelSD',ones(2)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeRelSD',-1),'matRad:Error');

function instanceTest_absRangeUncertainty(model)
    newValue = 2;
    model.rangeAbsSD = newValue;
    getValue = model.rangeAbsSD;
    assertEqual(newValue,getValue);
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeAbsSD','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeAbsSD',ones(2)),'matRad:Error'); 
    assertExceptionThrown(@() assignmentTestHelper(model,'rangeAbsSD',-1),'matRad:Error');

function instanceTest_shiftUncertainty(model)
    newValue = [1 1 1];
    model.shiftSD = newValue;
    getValue = model.shiftSD;
    assertEqual(newValue,getValue);
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD',ones(3,3)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD',ones(3,1)),'matRad:Error');    
    assertExceptionThrown(@() assignmentTestHelper(model,'shiftSD',[-1 2 2]),'matRad:Error');

function instanceTest_wcSigma(model)
    newValue = 5;
    model.wcSigma = newValue;
    getValue = model.wcSigma;
    assertEqual(newValue,getValue); 
    assertExceptionThrown(@() assignmentTestHelper(model,'wcSigma','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'wcSigma',-1),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'wcSigma',ones(2)),'matRad:Error');

function instanceTest_ctScenProb(model)
    newValue = rand(model.numOfCtScen,1);
    newValue = newValue ./ sum(newValue);
    model.ctScenProb(:,2) = newValue;
    getValue = model.ctScenProb(:,2);
    assertEqual(newValue,getValue); 
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb','a'),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb',ones(1,5)),'matRad:Error');
    assertExceptionThrown(@() assignmentTestHelper(model,'ctScenProb',-1*ones(model.numOfCtScen,1)),'matRad:Error');

function instanceTest_TYPE(model)
    %assertWarning(@() model.TYPE,'matRad:Deprecated');
    assertEqual(model.TYPE,model.name);

function instanceTest_wcFactor(model)
    %assertWarning(@() model.wcFactor,'matRad:Deprecated');
    assertEqual(model.TYPE,model.name);

