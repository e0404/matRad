function test_suite = test_interp3

test_functions = localfunctions();

initTestSuite;

function [xi, yi, zi, cube] = helper_cube()
xi = 0:2:6;
yi = 0:2:8;
zi = 0:2:4;
[X, Y, Z] = meshgrid(xi, yi, zi);
cube = X + 10 * Y + 100 * Z;

function test_defaultExtrapolationIsNaN
[xi, yi, zi, cube] = helper_cube();
res = matRad_interp3(xi, yi, zi, cube, [-1 1], [0 1]', [0 2]);
assertTrue(isnan(res(1, 1, 1)));
assertFalse(isnan(res(1, 2, 1)));

function test_numericExtrapolationValue
[xi, yi, zi, cube] = helper_cube();
res = matRad_interp3(xi, yi, zi, cube, [-1 1], [0 1]', [0 2], 'linear', -5);
assertEqual(res(1, 1, 1), -5);

function test_nearestExtrapolationReplicatesEdges
[xi, yi, zi, cube] = helper_cube();
xq = [-1 0 1 7];
yq = [-3 5 9]';
zq = [-1 5];

for mode = {'linear', 'nearest'}
    res = matRad_interp3(xi, yi, zi, cube, xq, yq, zq, mode{1}, 'nearest');
    assertTrue(all(isfinite(res(:))));

    % compare to the clamped query points
    ref = matRad_interp3(xi, yi, zi, cube, min(max(xq, 0), 6), min(max(yq, 0), 8), min(max(zq, 0), 4), mode{1});
    assertEqual(res, ref);
end

% corner values are those of the closest corner voxels
res = matRad_interp3(xi, yi, zi, cube, [-1 7], [-3 9]', [-1 5], 'linear', 'nearest');
assertEqual(res(1, 1, 1), cube(1, 1, 1));
assertEqual(res(end, end, end), cube(end, end, end));

function test_unknownExtrapolationModeThrows
[xi, yi, zi, cube] = helper_cube();
assertExceptionThrown(@() matRad_interp3(xi, yi, zi, cube, 1, 1, 1, 'linear', 'foo'), 'matRad:Error');
