function test_suite = test_collapsedConeTransport

test_functions = localfunctions();

initTestSuite;

function dirs = helper_randomDirections(n)
dirs = [0.3 -0.7 0.5; -0.8 0.1 0.6; 0.2 0.9 -0.4; 1 0 0; 0 0 -1; 0.577 0.577 0.577]';
dirs = repmat(dirs, 1, ceil(n / size(dirs, 2)));
dirs = dirs(:, 1:n) + 0.05 * reshape(sin(1:3 * n), 3, n);
dirs = dirs ./ repmat(sqrt(sum(dirs.^2, 1)), 3, 1);

function [terma, rho, voxelSize, dirs, wA, aA, wB, aB] = helper_heterogeneousInput()
dims = [13 17 11];
[i, j, k] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
terma = exp(-0.1 * i) .* (1 + 0.5 * sin(j)) .* (1 + 0.2 * cos(k));
rho = 0.6 + 0.5 * sin(0.7 * i + 0.3 * j) .* cos(0.5 * k);
rho(5:7, 6:9, 4:6) = 0;
voxelSize = [2 3 2.5];
dirs = helper_randomDirections(16);
wA = linspace(0.01, 0.05, 16);
aA = linspace(0.2, 1.5, 16);
wB = linspace(0.03, 0.01, 16);
aB = linspace(0.002, 0.02, 16);

function test_referenceDoseFiniteAndNonNegative
[terma, rho, voxelSize, dirs, wA, aA, wB, aB] = helper_heterogeneousInput();
dose = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, wA, aA, wB, aB);
assertEqual(size(dose), size(terma));
assertTrue(all(isfinite(dose(:))));
assertTrue(all(dose(:) >= 0));
assertTrue(any(dose(:) > 0));

function test_mexMatchesReference
if exist('matRad_collapsedConeMex', 'file') ~= 3
    moxunit_throw_test_skipped_exception('matRad_collapsedConeMex not compiled');
end
[terma, rho, voxelSize, dirs, wA, aA, wB, aB] = helper_heterogeneousInput();
doseRef = matRad_collapsedConeTransport(terma, rho, voxelSize, dirs, wA, aA, wB, aB);
doseMex = matRad_collapsedConeMex(terma, rho, voxelSize, dirs, wA, aA, wB, aB);
relTol = 1e-12;
assertElementsAlmostEqual(doseMex, doseRef, 'absolute', relTol * max(doseRef(:)));

function test_energyConservationShortRange
dims = [41 41 41];
terma = zeros(dims);
terma(16:26, 16:26, 16:26) = 1 + rand(11, 11, 11);
rho = ones(dims);
dirs = helper_randomDirections(12);
wA = 0.05 * ones(1, 12);
wB = 0.03 * ones(1, 12);
aA = 2 * ones(1, 12);
aB = 1 * ones(1, 12);
dose = matRad_collapsedConeTransport(terma, rho, [2 2 2], dirs, wA, aA, wB, aB);
assertElementsAlmostEqual(sum(dose(:)), sum(terma(:)) * sum(wA + wB), 'relative', 1e-6);

function test_equilibriumEachVoxelVisitedOnce
dims = [30 30 30];
terma = ones(dims);
rho = ones(dims);
dirs = helper_randomDirections(6);
for d = 1:size(dirs, 2)
    dose = matRad_collapsedConeTransport(terma, rho, [2 2 2], dirs(:, d), 0.1, 3, 0, 1);
    assertTrue(max(dose(:)) <= 0.1 * (1 + 1e-12));
    assertElementsAlmostEqual(median(dose(:)), 0.1, 'relative', 1e-6);
end

function test_densityScaling
dims = [60 5 5];
terma = zeros(dims);
terma(5, 3, 3) = 1;
dirs = [1; 0; 0];
doseWater = matRad_collapsedConeTransport(terma, ones(dims), [1 1 1], dirs, 0, 1, 1, 0.1);
doseHalf = matRad_collapsedConeTransport(terma, 0.5 * ones(dims), [2 1 1], dirs, 0, 1, 1, 0.1);
assertElementsAlmostEqual(doseHalf(:, 3, 3), doseWater(:, 3, 3), 'relative', 1e-12);

function test_axisDirectionMirrorSymmetry
dims = [25 21 19];
terma = zeros(dims);
terma(13, 11, 10) = 1;
rho = ones(dims);
for ax = 1:3
    axisDir = zeros(3, 1);
    axisDir(ax) = 1;
    doseForward = matRad_collapsedConeTransport(terma, rho, [2 3 2.5], axisDir, 0.1, 0.3, 0.05, 0.05);
    doseBackward = matRad_collapsedConeTransport(terma, rho, [2 3 2.5], -axisDir, 0.1, 0.3, 0.05, 0.05);
    assertElementsAlmostEqual(doseBackward, flip(doseForward, ax), 'absolute', 1e-15);
end

function test_fitCollapsedConeKernelSynthetic
edk.angleEdges = 0:3.75:180;
edk.radiusEdges = 10 * [0 0.05 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.8 1 1.5 2 3 4 5 6 8 10 15 20 30 40 50 60];
r1 = edk.radiusEdges(1:end - 1);
r2 = edk.radiusEdges(2:end);
theta = 0.5 * (edk.angleEdges(1:end - 1) + edk.angleEdges(2:end));
coneA = 0.02 * cosd(theta / 2).^4;
coneB = 0.01 * ones(size(theta));
edk.primary = coneA' * (exp(-0.3 * r1) - exp(-0.3 * r2));
edk.scatter = coneB' * (exp(-0.02 * r1) - exp(-0.02 * r2));

kernel = matRad_fitCollapsedConeKernel(edk, 12);
assertEqual(kernel.numPolarAngles, 12);
assertElementsAlmostEqual(kernel.polarEdges, 0:15:180);
assertElementsAlmostEqual(sum(kernel.A), sum(edk.primary(:)), 'relative', 1e-12);
assertElementsAlmostEqual(sum(kernel.B), sum(edk.scatter(:)), 'relative', 1e-12);
assertElementsAlmostEqual(kernel.a, 0.3 * ones(1, 12), 'relative', 1e-4);
assertElementsAlmostEqual(kernel.b, 0.02 * ones(1, 12), 'relative', 1e-3);
assertTrue(all(kernel.residual.primaryCum < 1e-4));

function test_fitCollapsedConeKernelInvalidBins
edk.angleEdges = 0:3.75:180;
edk.radiusEdges = [0 1 2];
edk.primary = ones(48, 2);
edk.scatter = ones(48, 2);
assertExceptionThrown(@() matRad_fitCollapsedConeKernel(edk, 7));
