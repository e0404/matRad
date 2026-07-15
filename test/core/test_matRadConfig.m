function test_suite = test_matRadConfig
%The output should always be test_suite, and the function name the same as
%your file name

test_functions=localfunctions();

initTestSuite;

function test_defaultsContainOxygen
    matRad_cfg = MatRad_Config.instance();

    % re-apply the standard defaults (this is what defines the oxygen
    % machine/bioModel defaults) and verify oxygen is present
    matRad_cfg.setDefaultProperties();

    assertTrue(isfield(matRad_cfg.defaults.machine, 'oxygen'));
    assertEqual(matRad_cfg.defaults.machine.oxygen, 'Generic');

    assertTrue(isfield(matRad_cfg.defaults.bioModel, 'oxygen'));
    assertEqual(matRad_cfg.defaults.bioModel.oxygen, 'LEM');

    % restore the testing defaults so subsequent tests are unaffected
    matRad_cfg.setDefaultPropertiesForTesting();

function test_logging_keepLogAndLevels
    matRad_cfg = MatRad_Config.instance();
    origKeepLog  = matRad_cfg.keepLog;
    origLogLevel = matRad_cfg.logLevel;

    % keeping the log stores messages in memory
    matRad_cfg.keepLog = true;
    matRad_cfg.logLevel = 5;
    matRad_cfg.dispInfo('test info message %d\n', 1);
    matRad_cfg.dispDebug('test debug message\n');
    matRad_cfg.dispDeprecationWarning('test deprecation\n');
    matRad_cfg.dispWarning('test warning\n');
    assertFalse(isempty(matRad_cfg.messageLog));

    % writing the kept log to a file
    logFile = [tempname '.log'];
    matRad_cfg.writeLogToFile(logFile);
    assertTrue(exist(logFile, 'file') == 2);
    delete(logFile);

    matRad_cfg.keepLog  = origKeepLog;
    matRad_cfg.logLevel = origLogLevel;

function test_logging_invalidTypeAndLevel
    matRad_cfg = MatRad_Config.instance();

    % an unknown log type raises an error
    assertExceptionThrown(@() matRad_cfg.displayToConsole('notAType', 'msg\n'));

    % an out-of-range log level is rejected by the set method
    origLogLevel = matRad_cfg.logLevel;
    assertExceptionThrown(@() setLogLevel(matRad_cfg, 99), 'matRad:Error');
    matRad_cfg.logLevel = origLogLevel;

function setLogLevel(cfg, value)
    cfg.logLevel = value;

function test_getDefaultProperties
    matRad_cfg = MatRad_Config.instance();

    % standard fields are populated from the defaults when missing
    pln = struct('radiationMode', 'protons');
    pln = matRad_cfg.getDefaultProperties(pln, 'propStf');
    assertTrue(isfield(pln, 'propStf'));

    % existing fields are merged with the defaults
    pln = matRad_cfg.getDefaultProperties(pln, {'propDoseCalc', 'propOpt'});
    assertTrue(isfield(pln, 'propDoseCalc'));
    assertTrue(isfield(pln, 'propOpt'));
