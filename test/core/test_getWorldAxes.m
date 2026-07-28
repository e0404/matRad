function test_suite = test_getWorldAxes

test_functions = localfunctions();

initTestSuite;

function test_getWorldAxesBasic
assertExceptionThrown(@()matRad_getWorldAxes());

ct = helper_getTestCt();
expectedX =  [-5 -4 -3 -2 -1 0 1 2 3 4];
expectedY = [-10 -8 -6 -4 -2 0 2 4 6 8];
expectedZ = [-15 -12 -9 -6 -3 0 3 6 9 12];
resultCt =  matRad_getWorldAxes(ct);
assertEqual(resultCt.x, expectedX);
assertEqual(resultCt.y, expectedY);
assertEqual(resultCt.z, expectedZ);

function test_getWorldAxesCompletesPartiallyPopulatedAxes
% Regression: the presence check reduced a 1x3 isfield result with an
% implicit all(), so it only detected the case where *every* axis was
% missing. A struct carrying just some of x/y/z fell through to the
% emptiness checks and errored on the first absent field.
ct = helper_getTestCt();
complete = matRad_getWorldAxes(ct);

partials = {struct('x', complete.x), ...
            struct('y', complete.y, 'z', complete.z), ...
            struct('x', complete.x, 'z', complete.z)};

for i = 1:numel(partials)
    gridStruct = ct;
    names = fieldnames(partials{i});
    for k = 1:numel(names)
        gridStruct.(names{k}) = partials{i}.(names{k});
    end

    result = matRad_getWorldAxes(gridStruct);
    assertEqual(result.x, complete.x);
    assertEqual(result.y, complete.y);
    assertEqual(result.z, complete.z);
end

function test_getWorldAxesRecomputesEmptyAxes
ct = helper_getTestCt();
complete = matRad_getWorldAxes(ct);

gridStruct = complete;
gridStruct.y = [];
result = matRad_getWorldAxes(gridStruct);

assertEqual(result.x, complete.x);
assertEqual(result.y, complete.y);
assertEqual(result.z, complete.z);

function test_getWorldAxesLeavesCompleteAxesUntouched
ct = helper_getTestCt();
complete = matRad_getWorldAxes(ct);

% populated axes are returned as-is, even when they disagree with
% resolution/cubeDim (callers may supply a shifted grid deliberately)
shifted = complete;
shifted.x = complete.x + 100;
result = matRad_getWorldAxes(shifted);

assertEqual(result.x, complete.x + 100);

function ct = helper_getTestCt()
ct.resolution.x = 1;
ct.resolution.y = 2;
ct.resolution.z = 3;
ct.cubeDim = [10 10 10];
