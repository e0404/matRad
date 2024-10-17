function test_suite = test_optimizerFmincon 
    
test_functions=localfunctions();

initTestSuite;

function test_optimizer_fmincon_construct

    if moxunit_util_platform_is_octave
        moxunit_throw_test_skipped_exception('fmincon not available for Octave!');
    end

    if ~license('test', 'optimization_toolbox') || license('checkout', 'optimization_toolbox') == 0
        moxunit_throw_test_skipped_exception('Optimization Toolbox containing fmincon not available!');
    end
    
    opti = matRad_OptimizerFmincon();
    assertTrue(isobject(opti));
    assertTrue(isa(opti, 'matRad_OptimizerFmincon'));
    
function test_optimizer_fmincon_available

    if moxunit_util_platform_is_octave
        moxunit_throw_test_skipped_exception('fmincon not available for Octave!');
    end

    if ~license('test', 'optimization_toolbox') || license('checkout', 'optimization_toolbox') == 0
        moxunit_throw_test_skipped_exception('Optimization Toolbox containing fmincon not available!');
    end

    opti = matRad_OptimizerFmincon();
    assertTrue(opti.IsAvailable());
    assertTrue(matRad_OptimizerFmincon.IsAvailable()); %Check static
    
function test_optimizer_fmincon_getStatus
    
        if moxunit_util_platform_is_octave
            moxunit_throw_test_skipped_exception('fmincon not available for Octave!');
        end
    
        if ~license('test', 'optimization_toolbox') || license('checkout', 'optimization_toolbox') == 0
            moxunit_throw_test_skipped_exception('Optimization Toolbox containing fmincon not available!');
        end
    
        opti = matRad_OptimizerFmincon();
        [statusmsg, statusflag] = opti.GetStatus();
        assertEqual(statusmsg, 'No Last Optimizer Status Available!');
        assertEqual(statusflag, -1);
        
        % TODO: test other status

% TODO: test optimize function