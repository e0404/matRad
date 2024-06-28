function test_suite = test_interp1

test_functions=localfunctions();

initTestSuite;

% Test basic values
function test_matrad_interp1_values
    %R = realmax; % For R = realmax, the test may often fail because of Inf in matRad_interp1: should first divide and then multiply.
    R = 10^100;

    % Pick a vector x and sort
    % First element of x is in [-R, R]
    x1el1 = (2*rand - 1)*R;
    if x1el1 < 0
        x1el2 = x1el1 + rand*(R);
    else
        x1el2 = x1el1 + rand*(R - x1el1);
    end
    if x1el1<=x1el2
        x1 = [x1el1; x1el2];
    else
        x1 = [x1el2; x1el1];
    end

    % Pick a sorted vector y
    y1el1 = (2*rand - 1)*R;
    if y1el1 < 0
        y1el2 = y1el1 + rand*(R);
    else
        y1el2 = y1el1 + rand*(R - y1el1);
    end
    y1 = [y1el1; y1el2];

    % Pick a single value intermediate in x
    x2 = rand*(x1(2)-x1(1)) + x1(1);

    % Verify interpolation works correctly
    y2 = matRad_interp1(x1, y1, x2);
    expectedy2 = y1(1) + (x2 - x1(1))*((y1(2) - y1(1))/(x1(2) - x1(1)));
    assertTrue(~isnan(y2));
    assertElementsAlmostEqual(y2, expectedy2);

    % Flip y vector and check interpolation again
    y1 = flip(y1);
    y2 = matRad_interp1(x1, y1, x2);
    expectedy2 = y1(1) + (x2 - x1(1))*(y1(2) - y1(1))/(x1(2) - x1(1));
    assertTrue(~isnan(y2));
    assertElementsAlmostEqual(y2, expectedy2);

    % Pick a x value exceeding upper boundaries
    % If x2 is outside boundaries, we expect a NaN
    x2 = x1(2) + rand*(R - x1(2));
    if x1(2) == R
        assertTrue(~isnan(matRad_interp1(x1, y1, x2)));
        assertEqual(x2, x1(2));
    else
        assertTrue(isnan(matRad_interp1(x1, y1, x2)));
    end
    % Pick a x value exceeding lower boundaries
    % If x2 is outside boundaries, we expect a NaN
    x2 = x1(1) + rand*(-R - x1(1));
    if x1(1) == -R
        assertTrue(~isnan(matRad_interp1(x1, y1, x2)));
        assertEqual(x2, x1(2));
    else
        assertTrue(isnan(matRad_interp1(x1, y1, x2)));
    end

    
% Test Extrapolation Methods
function test_matRad_interp1_extrapolation
    R = 10^100;                 % Choose a maximum scale
    realExtrap = R*(2*rand-1);  % Choose a random value for Real Extrapolation
   
    % Pick a vector x and sort
    % First element of x is in [-R, R]
    x1el1 = (2*rand - 1)*R;
    if x1el1 < 0
        x1el2 = x1el1 + rand*(R);
    else
        x1el2 = x1el1 + rand*(R - x1el1);
    end
    if x1el1<=x1el2
        x1 = [x1el1; x1el2];
    else
        x1 = [x1el2; x1el1];
    end

    % Pick a sorted vector y
    y1el1 = (2*rand - 1)*R;
    if y1el1 < 0
        y1el2 = y1el1 + rand*(R);
    else
        y1el2 = y1el1 + rand*(R - y1el1);
    end
    y1 = [y1el1; y1el2]; 

    % Pick 2 values exceeding Lower and Upper boundaries
    % Pick 1 value in the boundaries
    xLow = x1(1) + rand*(-R - x1(1));
    xUp = x1(2) + rand*(R - x1(2));
    x2 = rand*(x1(2)-x1(1)) + x1(1);
    x2 = [xLow; x2; xUp];
    
    % Index vector for out-of-boundaries values
    outIdx = find(x2<x1(1) | x2>x1(2));
    % [1] Real extrapolation: we expect a constant value for
    % out-of-bondaries interpolation
    y2 = matRad_interp1(x1, y1, x2, realExtrap);
    assertElementsAlmostEqual(y2(outIdx), realExtrap.*ones(size(y2(outIdx))));
    % [2] NaN & 'none': we expect NaNs for both of them in this case
    % If numel(x2) == 1, we expect an error for 'none'----> Implement this!
    y2 = matRad_interp1(x1, y1, x2, NaN);
    y3 = matRad_interp1(x1, y1, x2, 'none');
    assertTrue( sum(isnan(y2(outIdx)))==length(outIdx) );
    assertTrue( sum(isnan(y3(outIdx)))==length(outIdx) );
    assertElementsAlmostEqual(y2, y3);
    % [3] linear & extrap: in this case we expect they work the same
    % If numel(x2) == 1, we expect an error for 'linear'----> Implement this!
    y2 = matRad_interp1(x1, y1, x2, 'extrap');
    y3 = matRad_interp1(x1, y1, x2, 'linear'); 
    assertElementsAlmostEqual(y2, y3);

function test_matRad_interp1_errors
    R = 10^100;

    % Repetition Errors
    % Pick x1 with repetitions and y1 sorted
    x1 = (2*rand - 1)*R.*[1;1];
    y1el1 = (2*rand - 1)*R;
    if y1el1 < 0
        y1el2 = y1el1 + rand*(R);
    else
        y1el2 = y1el1 + rand*(R - y1el1);
    end
    y1 = [y1el1; y1el2];
    % [1] x2 is a vector: we expect error
    % First value is repeted value in x1, second value is different; 
    x2 = [x1(1); rand.*((R-x1(1))/10)+x1(1)];
    assertExceptionThrown(@() matRad_interp1(x1, y1, x2), 'MATLAB:griddedInterpolant:NonUniqueCompVecsPtsErrId');

    % [2] x2 is a scalar, the result is NaN
    y2 = matRad_interp1(x1, y1, x2(1));
    assertTrue(isnan(y2));
    y2 = matRad_interp1(x1, y1, x2(2));
    assertTrue(isnan(y2));

    % Extrapolation Errors
    % [1] Single query point, case 'none'
    x1 = [1; 10000];
    y1 = [7.3; 2.4];
    assertExceptionThrown(@() matRad_interp1(x1, y1, 0.5, 'none'), '');

    % [2] Single query point, case 'linear'
    assertExceptionThrown(@() matRad_interp1(x1, y1, 0.5, 'linear'), '');

%{
function test_matRad_interp1_multiple1D
    R = 10^100;
    x1 = (2.*rand(1000, 1)- 1).*R ;
    x1 = unique(x1);
    x1 = sort(x1);
    y1 = (2.*rand(size(x1))- 1).*R;
%}


