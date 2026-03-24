function test_suite = test_rayTracer

test_functions = localfunctions();

initTestSuite;

function test_siddonRayTracer
% test function with dummy numerical example

cubes{1} = ones([2, 2, 2]);
cubes{2} = cubes{1};
cubes{2}(:, :, 2) = [2, 2; 2, 2];

resolution.x = 1;
resolution.y = 1;
resolution.z = 1;

grid.resolution = resolution;
grid.dimensions = size(cubes{1});

isocenter   = [0, 0, 0];
sourcePoint = [-1, -1, -2];
targetPoint = [0,  0, 2];

% Now we will have voxels centers at -1 0 (default matRad world
% coordinates, meaning the isocenter will point into the [2 2 2] voxel)
% A ray starting at -1 -1 -2 and ending at 0 0 1 sees the first plane at
% z = -1.5 and the last plae at z=0.5. When intersecting the first
% plane, it will have passed 1/8 of its length.

rt = matRad_RayTracerSiddon(cubes, grid);
[alphas, l, rho, d12, ix] = rt.traceRay(isocenter, sourcePoint, targetPoint);

% test Output types
assertTrue(isvector(alphas));
assertTrue(isvector(l));
assertTrue(iscell(rho));
assertTrue(isfloat(d12));
assertTrue(isvector(ix));

% test numerical Output
grid = matRad_getWorldAxes(grid);
rayVec = targetPoint - sourcePoint;
rayLength = norm(rayVec);
% the ray will intersect z at the coordinates of the z voxel boundaries

entryPoints = [sourcePoint + rayVec * 1 / 8
               sourcePoint + rayVec * 3 / 8
               sourcePoint + rayVec * 5 / 8];
entryPoints = entryPoints - sourcePoint;
alphasNum   = sqrt(sum(entryPoints.^2, 2)) ./ rayLength;
lNum        = [rayLength / 4, rayLength / 4];
rhoNum{1}   = [1, 1];
rhoNum{2}   = [1, 2];
d12Num      = rayLength;
ixNum       = [1, 8];

assertElementsAlmostEqual(alphasNum', alphas);
assertElementsAlmostEqual(lNum, l);
assertElementsAlmostEqual(rhoNum{1}, rho{1});
assertElementsAlmostEqual(rhoNum{2}, rho{2});
assertElementsAlmostEqual(d12Num, d12);
assertEqual(ixNum, ix);

% test the old deprecated function with dummy numerical example
% It expects cube coords for the isocenter

isocenterCube = matRad_world2cubeCoords(isocenter, grid);
[alphasOld, lOld, rhoOld, d12Old, ixOld] = matRad_siddonRayTracer(isocenterCube, ...
                                                                  resolution, ...
                                                                  sourcePoint, ...
                                                                  targetPoint, ...
                                                                  cubes);

assertElementsAlmostEqual(alphas, alphasOld);
assertElementsAlmostEqual(l, lOld);
assertElementsAlmostEqual(rho{1}, rhoOld{1});
assertElementsAlmostEqual(rho{2}, rhoOld{2});
assertElementsAlmostEqual(d12, d12Old);
assertEqual(ixNum, ixOld);

function test_rayDoesNotHitCT

cubes{1} = ones([2, 2, 2]);

resolution.x = 1;
resolution.y = 1;
resolution.z = 1;

grid.resolution = resolution;
grid.dimensions = size(cubes{1});

isocenter   = [-2, -2, -2];
sourcePoint = [2.5, 2.5 -4];
targetPoint = [10, 10, 6];

rt = matRad_RayTracerSiddon(cubes, grid);
[alphas, l, rho, d12, ix] = rt.traceRay(isocenter, sourcePoint, targetPoint);

% test numerical Output
d12Num      = norm(sourcePoint - targetPoint);

assertTrue(isempty(alphas));
assertTrue(isempty(l));
assertTrue(isempty(ix));
assertTrue(isempty(rho{1}));
assertElementsAlmostEqual(d12Num, d12);

% deprecated call using cube coordinates
isocenterCube = matRad_world2cubeCoords(isocenter, grid);
[alphasOld, lOld, rhoOld, d12Old, ixOld] = matRad_siddonRayTracer(isocenterCube, ...
                                                                  resolution, ...
                                                                  sourcePoint, ...
                                                                  targetPoint, ...
                                                                  cubes);

assertElementsAlmostEqual(alphas, alphasOld);
assertElementsAlmostEqual(l, lOld);
assertElementsAlmostEqual(rho{1}, rhoOld{1});
assertElementsAlmostEqual(d12, d12Old);
assertEqual(ix, ixOld);

function test_rayHitsAtBoundary

cubes{1} = ones([2, 2, 2]);
cubes{2} = cubes{1};
cubes{2}(:, :, 2) = [2, 2; 2, 2];

resolution.x = 1;
resolution.y = 1;
resolution.z = 1;

isocenter   = [-2, -2, -2];
sourcePoint = [2.5, 2.5, -4];
targetPoint = [2.5, 2.5, 6];

grid.resolution = resolution;
grid.dimensions = size(cubes{1});

rt = matRad_RayTracerSiddon(cubes, grid);
[alphas, l, rho, d12, ix] = rt.traceRay(isocenter, sourcePoint, targetPoint);

% test Output types
assertTrue(isvector(alphas));
assertTrue(isvector(l));
assertTrue(iscell(rho));
assertTrue(isfloat(d12));
assertTrue(isvector(ix));

% test numerical Output
entryPoints = [2.5, 2.5, 0.5
               2.5, 2.5, 1.5
               2.5, 2.5, 2.5];
entryPoints = entryPoints - sourcePoint;
alphasNum   = sqrt(sum(entryPoints.^2, 2)) ./ 10;
lNum        = [1, 1];
rhoNum{1}   = [1, 1];
rhoNum{2}   = [1, 2];
d12Num      = 10;
ixNum       = [4, 8];

assertElementsAlmostEqual(alphasNum', alphas);
assertElementsAlmostEqual(lNum, l);
assertElementsAlmostEqual(rhoNum{1}, rho{1});
assertElementsAlmostEqual(rhoNum{2}, rho{2});
assertElementsAlmostEqual(d12Num, d12);
assertEqual(ixNum, ix);

% deprecated call using cube coordinates
isocenterCube = matRad_world2cubeCoords(isocenter, grid);
[alphasOld, lOld, rhoOld, d12Old, ixOld] = matRad_siddonRayTracer(isocenterCube, ...
                                                                  resolution, ...
                                                                  sourcePoint, ...
                                                                  targetPoint, ...
                                                                  cubes);

assertElementsAlmostEqual(alphas, alphasOld);
assertElementsAlmostEqual(l, lOld);
assertElementsAlmostEqual(rho{1}, rhoOld{1});
assertElementsAlmostEqual(rho{2}, rhoOld{2});
assertElementsAlmostEqual(d12, d12Old);
assertEqual(ix, ixOld);

function test_rayHitsAtCorner

cubes{1} = ones([2, 2, 2]);
cubes{2} = cubes{1};
cubes{2}(:, :, 2) = [2, 2; 2, 2];

resolution.x = 1;
resolution.y = 1;
resolution.z = 1;

isocenter   = [-2, -2, -2];
sourcePoint = [1.5, 1.5, -4];
targetPoint = [3.5, 3.5, 5];

grid.resolution = resolution;
grid.dimensions = size(cubes{1});

rt = matRad_RayTracerSiddon(cubes, grid);
[alphas, l, rho, d12, ix] = rt.traceRay(isocenter, sourcePoint, targetPoint);

% test numerical Output
alphasNum   = 0.5;
ixNum       = 1:0;
rhoNum{1}   = cubes{1}(ixNum);
rhoNum{2}   = cubes{2}(ixNum);
d12Num      = norm(sourcePoint - targetPoint);

assertElementsAlmostEqual(alphasNum', alphas);
assertTrue(isempty(l));
assertEqual(size(l), [1 0]);
assertElementsAlmostEqual(rhoNum{1}, rho{1});
assertElementsAlmostEqual(rhoNum{2}, rho{2});
assertElementsAlmostEqual(d12Num, d12);
assertEqual(ixNum, ix);

% deprecated call using cube coordinates
isocenterCube = matRad_world2cubeCoords(isocenter, grid);
[alphasOld, lOld, rhoOld, d12Old, ixOld] = matRad_siddonRayTracer(isocenterCube, ...
                                                                  resolution, ...
                                                                  sourcePoint, ...
                                                                  targetPoint, ...
                                                                  cubes);

assertElementsAlmostEqual(alphas, alphasOld);
assertElementsAlmostEqual(l, lOld);
assertElementsAlmostEqual(rho{1}, rhoOld{1});
assertElementsAlmostEqual(rho{2}, rhoOld{2});
assertElementsAlmostEqual(d12, d12Old);
assertElementsAlmostEqual(ix, ixOld);

function test_vectorizedVsLoop

testData = load('photons_testData.mat');
targetPoints = vertcat(testData.stf(1).ray.targetPoint);
sourcePoint = testData.stf(1).sourcePoint;
isocenter = testData.stf(1).isoCenter;

rt = matRad_RayTracerSiddon(testData.ct.cube, testData.ct);

rt.vectorized = true;
[alphas, l, rho, d12, ix] = rt.traceRays(isocenter, sourcePoint, targetPoints);

rt.vectorized = false;
[alphasLoop, lLoop, rhoLoop, d12Loop, ixLoop] = rt.traceRays(isocenter, sourcePoint, targetPoints);

assertElementsAlmostEqual(alphas, alphasLoop);
assertElementsAlmostEqual(l, lLoop);
assertElementsAlmostEqual(d12, d12Loop);
assertElementsAlmostEqual(ix, ixLoop);
assertElementsAlmostEqual(rho{1}, rhoLoop{1});

function test_singleVsDouble
testData = load('photons_testData.mat');
targetPoints = vertcat(testData.stf(1).ray.targetPoint);
sourcePoint = testData.stf(1).sourcePoint;
isocenter = testData.stf(1).isoCenter;

rt = matRad_RayTracerSiddon(testData.ct.cube, testData.ct);

% Test force double
rt.forcePrecision = 'double';
[alphas, l, rho, d12, ix] = rt.traceRays(isocenter, sourcePoint, targetPoints);

assertTrue(isa(alphas, 'double'));
assertTrue(isa(l, 'double'));
assertTrue(isa(d12, 'double'));
assertTrue(isa(ix, 'double'));
assertTrue(isa(rho{1}, 'double'));

% test force single
rt.forcePrecision = 'single';
[alphasSingle, lSingle, rhoSingle, d12Single, ixSingle] = rt.traceRays(isocenter, sourcePoint, targetPoints);

assertTrue(isa(alphasSingle, 'single'));
assertTrue(isa(lSingle, 'single'));
assertTrue(isa(d12Single, 'single'));
assertTrue(isa(ixSingle, 'double'));
assertTrue(isa(rhoSingle{1}, 'single'));

assertElementsAlmostEqual(alphasSingle, alphas);
assertElementsAlmostEqual(lSingle, l);
assertElementsAlmostEqual(d12Single, d12);
assertElementsAlmostEqual(ixSingle, ix);
assertElementsAlmostEqual(rhoSingle{1}, rho{1});
