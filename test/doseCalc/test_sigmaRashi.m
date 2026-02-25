function test_suite= test_sigmaRashi

test_functions=localfunctions();

initTestSuite;

function test_calcSigmaRashi

    baseDataEntry.range = 100;
    
    rangeShifter.ID = 1;
    rangeShifter.eqThickness = 1;
    rangeShifter.sourceRashiDistance = 9000;

    SSD = 10000;

    sigma = matRad_calcSigmaRashi(baseDataEntry, rangeShifter, SSD);

    assertElementsAlmostEqual(sigma,1.1,'relative',1e-2,1e-2);
