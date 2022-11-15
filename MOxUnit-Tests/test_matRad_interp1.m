function test_suite=test_matRad_interp1
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function testInterp1Extrapolation()
    xi = [100, 150, 200, 250, 300, 400, 500];
    yi = [2506.7, 2582.8, 2658.1, 2733.7, 2810.4, 2967.9, 3131.6]';
    x = 650;
    matRad_interp1_solution = matRad_interp1(xi, yi, x);
    interp1_solution = interp1(xi, yi, x);
    assertEqual(matRad_interp1_solution,interp1_solution)

function testInterp1Interpolation()
    xi = [100, 150, 200, 250, 300, 400, 500];
    yi = [2506.7, 2582.8, 2658.1, 2733.7, 2810.4, 2967.9, 3131.6]';
    x = 222;
    matRad_interp1_solution = matRad_interp1(xi, yi, x);
    interp1_solution = interp1(xi, yi, x);
    assertEqual(matRad_interp1_solution,interp1_solution)
    
function testFalseInput()
    xi = [100, 150, 200, 300, 400, 500]; % xi.length < yi.length => false input
    yi = [2506.7, 2582.8, 2658.1, 2733.7, 2810.4, 2967.9, 3131.6]';
    x = 2680.78;
    assertExceptionThrown(@()interp1(xi, yi, x));