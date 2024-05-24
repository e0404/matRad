function test_suite = test_scenarioModel

test_functions=localfunctions();

initTestSuite;

function test_scenarioAbstract
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@() matRad_ScenarioModel(),'');
    else
        assertExceptionThrown(@() matRad_ScenarioModel(),'MATLAB:class:abstract');
    end

