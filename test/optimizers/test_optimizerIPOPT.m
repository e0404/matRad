function test_suite = test_optimizerIPOPT 

test_functions=localfunctions();

initTestSuite;

function test_optimizer_ipopt_construct

    opti = matRad_OptimizerIPOPT();
    assertTrue(isobject(opti));
    assertTrue(isa(opti, 'matRad_OptimizerIPOPT'));
    
function test_optimizer_ipopt_available

    opti = matRad_OptimizerIPOPT();
    assertTrue(opti.IsAvailable());
    assertTrue(matRad_OptimizerIPOPT.IsAvailable()); %Check static
    
function test_optimizer_ipopt_getStatus

    opti = matRad_OptimizerIPOPT();
    [statusmsg, statusflag] = opti.GetStatus();
    assertEqual(statusmsg, 'No Last IPOPT Status Available!');
    assertEqual(statusflag, -1);
        
    % TODO: test other status

% TODO: test optimize function
