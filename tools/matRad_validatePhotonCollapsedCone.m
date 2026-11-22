function results = matRad_validatePhotonCollapsedCone(varargin)
% matRad_validatePhotonCollapsedCone compares the collapsed cone engine
% (CCC) with the SVD pencil beam engine (SVDPB) and, if available, the
% ompMC Monte Carlo engine
%
% Cases:
%   water   300 mm water cube, 10 x 10 cm field: central axis depth dose
%           and lateral profiles at 15, 50, 100 and 200 mm depth
%   slab10  water (0-50 mm) / lung (rED 0.25, 50-150 mm) / water slab
%           phantom, 10 x 10 cm field: central axis depth dose
%   slab3   same slab phantom, 3 x 3 cm field
% For the slab cases the central axis dose in the lung and at the distal
% lung/water interface is compared with ompMC. For all cases the roughness
% of the central axis dose (rms of the second difference divided by the
% mean) is reported per homogeneous region, excluding 10 mm around
% interfaces and the phantom boundaries, for the kernel tilts 'beam' and
% 'ray'.
%   fields  water cube, central axis dose at 100 mm depth for 3 x 3,
%           5 x 5, 10 x 10 and 20 x 20 cm fields, normalized to 10 x 10, and
%           the dose at 100 mm depth 30 mm outside the 10 x 10 field edge,
%           for CCC, SVDPB and ompMC
% The convergence of
% the CCC dose with respect to scatterCutOff is studied for a 20 x 20 cm
% field in water.
% Computation times are reported for all engines.
%
% call:
%   results = matRad_validatePhotonCollapsedCone()
%   results = matRad_validatePhotonCollapsedCone('Name', value, ...)
%
% input (name-value pairs):
%   resolution:     CT and dose grid resolution [mm]. Default: 4
%   scatterCutOffs: scatterCutOff values of the convergence study [mm].
%                   Default: [60 100 150 200]
%   runOmpMC:       run the ompMC engine if available. Default: true
%   numHistories:   total ompMC histories per case (photon splitting 20
%                   is applied by the ompMC engine). Default: 2e6
%                   ompMC is run as dose influence calculation and the dose
%                   is obtained as dij * w, since its forward calculation
%                   does not support more than one bixel per beam.
%   fieldSizeHistories: ompMC histories for the 3, 5, 10 and 20 cm fields
%                   of the field size study. Default: [5e5 5e5 2e6 8e6]
%   mcRelDoseThreshold: relative dose threshold of the ompMC beamlets
%                   (the matRad default of 0.005 removes the scatter tails
%                   of each beamlet). Default: 1e-5
%   averageWidth:   the central axis dose is averaged over
%                   (2 * averageWidth + 1)^2 voxels perpendicular to the
%                   beam for all engines, to reduce the Monte Carlo noise.
%                   Default: 1
%   plot:           show figures. Default: true
%
% output:
%   results:        struct with depth dose curves, profiles, roughness,
%                   convergence and timing data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser();
addParameter(p, 'resolution', 4);
addParameter(p, 'scatterCutOffs', [60 100 150 200]);
addParameter(p, 'runOmpMC', true);
addParameter(p, 'numHistories', 2e6);
addParameter(p, 'averageWidth', 1);
addParameter(p, 'mcRelDoseThreshold', 1e-5);
addParameter(p, 'fieldSizeHistories', [5e5 5e5 2e6 8e6]);
addParameter(p, 'plot', true);
parse(p, varargin{:});
opts = p.Results;

matRad_cfg = matRad_rc();
logLevel = matRad_cfg.logLevel;
matRad_cfg.logLevel = 2;
restoreLog = onCleanup(@() matRad_setLogLevel(matRad_cfg, logLevel));

lungDensity = 0.25;
results.water = matRad_validationCase(opts, [], 100, 20);
results.slab10 = matRad_validationCase(opts, [50 150 lungDensity], 100, 20);
results.slab3 = matRad_validationCase(opts, [50 150 lungDensity], 30, 10);
results.fields = matRad_validationFieldSizes(opts);
results.convergence = matRad_validationConvergence(opts);

matRad_cfg.logLevel = max(logLevel, 3);
matRad_printValidation(results);
if opts.plot
    matRad_plotValidation(results);
end

end

function matRad_setLogLevel(cfg, logLevel)
cfg.logLevel = logLevel;
end

function [ct, cst, pln, stf] = matRad_validationPhantom(opts, slab, fieldSize, bixelWidth)
res = opts.resolution;
n = round(300 / res);
builder = matRad_PhantomBuilder([n n n], [res res res], 1);
builder.addBoxTarget('Target', round(max(fieldSize + 20, 60) / res) * [1 1 1], 'HU', 0, ...
                     'objectives', struct(DoseObjectives.matRad_SquaredDeviation(800, 60)));
builder.addBoxOAR('Water', (n - 4) * [1 1 1], 'HU', 0, ...
                  'objectives', struct(DoseObjectives.matRad_SquaredOverdosing(10, 30)));
[ct, cst] = builder.getctcst();

if ~isempty(slab)
    rho = ones(ct.cubeDim);
    rho(ct.cubeHU{1} <= -1000) = 0.001;
    first = find(ct.cubeHU{1}(:, round(n / 2), round(n / 2)) > -1000, 1);
    lungRows = first + (round(slab(1) / res):round(slab(2) / res) - 1);
    rho(lungRows, 3:n - 2, 3:n - 2) = slab(3);
    hlut = matRad_loadHLUT(ct, 'photons');
    [rhoTable, ix] = unique(hlut(:, 2));
    ct.cubeHU = {reshape(interp1(rhoTable, hlut(ix, 1), rho(:), 'linear', 'extrap'), ct.cubeDim)};
end

pln.radiationMode = 'photons';
pln.machine = 'Generic';
pln.numOfFractions = 1;
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.bixelWidth = bixelWidth;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propDoseCalc.doseGrid.resolution = struct('x', res, 'y', res, 'z', res);
stf = matRad_generateStf(ct, cst, pln);

rayPos = reshape([stf.ray.rayPos_bev], 3, [])';
keep = all(abs(rayPos(:, [1 3])) <= (fieldSize - bixelWidth) / 2 + 1e-6, 2);
stf.ray = stf.ray(keep);
stf.numOfRays = nnz(keep);
stf.numOfBixelsPerRay = ones(1, stf.numOfRays);
stf.totalNumOfBixels = stf.numOfRays;
end

function out = matRad_validationCase(opts, slab, fieldSize, bixelWidth)
matRad_cfg = MatRad_Config.instance();
[ct, cst, pln, stf] = matRad_validationPhantom(opts, slab, fieldSize, bixelWidth);
out.fieldSize = fieldSize;
out.slab = slab;
w = ones(stf.totalNumOfBixels, 1);
center = round(ct.cubeDim / 2);
res = opts.resolution;

engines = matRad_validationEngines(opts);
out.engines = fieldnames(engines);
for e = 1:numel(out.engines)
    name = out.engines{e};
    plnE = pln;
    plnE.propDoseCalc = engines.(name);
    plnE.propDoseCalc.doseGrid = pln.propDoseCalc.doseGrid;
    matRad_cfg.dispInfo('Running %s ...\n', name);
    tStart = tic;
    resultGUI = matRad_validationForward(ct, cst, stf, plnE, w);
    out.time.(name) = toc(tStart);
    dose = resultGUI.physicalDose;
    isBody = ct.cubeHU{1}(:, center(2), center(3)) > -1000;
    first = find(isBody, 1);
    lat = center(2) + (-opts.averageWidth:opts.averageWidth);
    lon = center(3) + (-opts.averageWidth:opts.averageWidth);
    column = mean(reshape(dose(:, lat, lon), size(dose, 1), []), 2);
    if isfield(resultGUI, 'physicalDose_MCvar')
        variance = sum(reshape(resultGUI.physicalDose_MCvar(:, lat, lon), size(dose, 1), []), 2) / numel(lat)^2 / numel(lon)^2;
        out.relStd.(name) = sqrt(variance(first:end)) ./ column(first:end);
    end
    out.exitDepth = (find(isBody, 1, 'last') - first + 1) * res;
    out.depth = ((first:ct.cubeDim(1)) - first + 0.5)' * res;
    out.pdd.(name) = column(first:end);
    rowIx = sub2ind(ct.cubeDim, center(1) * ones(ct.cubeDim(2), 1), (1:ct.cubeDim(2))', center(3) * ones(ct.cubeDim(2), 1));
    worldCoords = matRad_cubeIndex2worldCoords(rowIx, ct);
    out.lateral = worldCoords(:, 1) - pln.propStf.isoCenter(1);
    for d = [15 50 100 200]
        row = first + round(d / res - 0.5);
        out.profiles.(name).(sprintf('d%d', d)) = squeeze(dose(row, :, center(3)))';
    end
end
out.roughness = matRad_validationRoughness(out, slab);
end

function resultGUI = matRad_validationForward(ct, cst, stf, pln, w)
if ~strcmp(pln.propDoseCalc.engine, 'ompMC')
    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);
    return
end
pln.propDoseCalc.numHistoriesPerBeamlet = ceil(pln.propDoseCalc.numHistoriesDirect / numel(w));
pln.propDoseCalc = rmfield(pln.propDoseCalc, 'numHistoriesDirect');
pln.propDoseCalc.outputMCvariance = true;
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);
dose = reshape(full(dij.physicalDose{1} * w), dij.doseGrid.dimensions);
variance = reshape(full(dij.physicalDose_MCvar{1} * w.^2), dij.doseGrid.dimensions);
resultGUI.physicalDose = matRad_interp3(dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z, dose, ...
                                        dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z, 'linear', 0);
resultGUI.physicalDose_MCvar = matRad_interp3(dij.doseGrid.x, dij.doseGrid.y', dij.doseGrid.z, variance, ...
                                              dij.ctGrid.x, dij.ctGrid.y', dij.ctGrid.z, 'linear', 0);
end

function engines = matRad_validationEngines(opts)
engines.SVDPB = struct('engine', 'SVDPB', 'kernelCutOff', Inf, 'geometricLateralCutOff', 180, 'enableDijSampling', false);
engines.CCC = struct('engine', 'CCC', 'enableDijSampling', false);
engines.CCCray = struct('engine', 'CCC', 'enableDijSampling', false, 'kernelTilt', 'ray');
if opts.runOmpMC
    plnMC = struct('radiationMode', 'photons', 'machine', 'Generic');
    if DoseEngines.matRad_PhotonOmpMCEngine.isAvailable(plnMC) && exist('omc_matrad', 'file') == 3
        engines.ompMC = struct('engine', 'ompMC', 'numHistoriesDirect', opts.numHistories, ...
                               'relativeDosimetricCutOff', 1 - opts.mcRelDoseThreshold);
    end
end
end

function roughness = matRad_validationRoughness(out, slab)
margin = 10;
if isempty(slab)
    edges = [0 out.exitDepth];
else
    edges = [0 slab(1) slab(2) out.exitDepth];
end
for e = 1:numel(out.engines)
    name = out.engines{e};
    r = zeros(1, numel(edges) - 1);
    for i = 1:numel(edges) - 1
        sel = out.depth >= edges(i) + margin & out.depth <= edges(i + 1) - margin;
        p = out.pdd.(name)(sel);
        r(i) = sqrt(mean(diff(p, 2).^2)) / mean(p);
    end
    roughness.(name) = r;
end
end

function fs = matRad_validationFieldSizes(opts)
matRad_cfg = MatRad_Config.instance();
fs.fieldSize = [30 50 100 200];
bixelWidths = [10 10 20 20];
outOfField = 80;
engines = matRad_validationEngines(opts);
engines = rmfield(engines, 'CCCray');
fs.engines = fieldnames(engines);
for f = 1:numel(fs.fieldSize)
    [ct, cst, pln, stf] = matRad_validationPhantom(opts, [], fs.fieldSize(f), bixelWidths(f));
    for e = 1:numel(fs.engines)
        name = fs.engines{e};
        plnE = pln;
        plnE.propDoseCalc = engines.(name);
        plnE.propDoseCalc.doseGrid = pln.propDoseCalc.doseGrid;
        if strcmp(name, 'ompMC')
            plnE.propDoseCalc.numHistoriesDirect = opts.fieldSizeHistories(f);
        end
        matRad_cfg.dispInfo('Field size %g mm: running %s ...\n', fs.fieldSize(f), name);
        tStart = tic;
        resultGUI = matRad_validationForward(ct, cst, stf, plnE, ones(stf.totalNumOfBixels, 1));
        fs.time.(name)(f) = toc(tStart);
        [fs.cax.(name)(f), fs.caxRelStd.(name)(f)] = matRad_validationPointDose(resultGUI, ct, opts, 100, 0);
        if fs.fieldSize(f) == 100
            [fs.out.(name), fs.outRelStd.(name)] = matRad_validationPointDose(resultGUI, ct, opts, 100, outOfField);
        end
    end
end
end

function [dose, relStd] = matRad_validationPointDose(resultGUI, ct, opts, depth, lateralOffset)
res = opts.resolution;
center = round(ct.cubeDim / 2);
first = find(ct.cubeHU{1}(:, center(2), center(3)) > -1000, 1);
lat = center(2) + round(lateralOffset / res) + (-opts.averageWidth:opts.averageWidth);
lon = center(3) + (-opts.averageWidth:opts.averageWidth);
numVox = numel(lat) * numel(lon);
rows = first:ct.cubeDim(1);
depths = (rows - first + 0.5)' * res;
column = mean(reshape(resultGUI.physicalDose(rows, lat, lon), numel(rows), []), 2);
dose = interp1(depths, column, depth);
relStd = NaN;
if isfield(resultGUI, 'physicalDose_MCvar')
    variance = sum(reshape(resultGUI.physicalDose_MCvar(rows, lat, lon), numel(rows), []), 2) / numVox^2;
    relStd = interp1(depths, sqrt(variance), depth) / dose;
end
end

function conv = matRad_validationConvergence(opts)
matRad_cfg = MatRad_Config.instance();
[ct, cst, pln, stf] = matRad_validationPhantom(opts, [], 200, 20);
w = ones(stf.totalNumOfBixels, 1);
center = round(ct.cubeDim / 2);
first = find(ct.cubeHU{1}(:, center(2), center(3)) > -1000, 1);
depth = ((first:ct.cubeDim(1)) - first + 0.5)' * opts.resolution;

conv.scatterCutOff = opts.scatterCutOffs;
conv.dose100 = zeros(size(conv.scatterCutOff));
conv.dose200 = zeros(size(conv.scatterCutOff));
conv.time = zeros(size(conv.scatterCutOff));
for i = 1:numel(conv.scatterCutOff)
    pln.propDoseCalc = struct('engine', 'CCC', 'enableDijSampling', false, 'scatterCutOff', conv.scatterCutOff(i), ...
                              'doseGrid', pln.propDoseCalc.doseGrid);
    matRad_cfg.dispInfo('Convergence: scatterCutOff = %g mm\n', conv.scatterCutOff(i));
    tStart = tic;
    resultGUI = matRad_calcDoseForward(ct, cst, stf, pln, w);
    conv.time(i) = toc(tStart);
    column = resultGUI.physicalDose(first:end, center(2), center(3));
    conv.dose100(i) = interp1(depth, column, 100);
    conv.dose200(i) = interp1(depth, column, 200);
end
end

function matRad_printValidation(results)
matRad_cfg = MatRad_Config.instance();
cases = {'water', 'slab10', 'slab3'};
for c = 1:numel(cases)
    out = results.(cases{c});
    matRad_cfg.dispInfo('\n=== %s: central axis dose ===\n%8s', cases{c}, 'd [mm]');
    matRad_cfg.dispInfo(' %9s', out.engines{:});
    matRad_cfg.dispInfo('\n');
    for d = [5 15 30 50 75 100 125 150 175 200 250]
        matRad_cfg.dispInfo('%8g', d);
        for e = 1:numel(out.engines)
            matRad_cfg.dispInfo(' %9.4f', interp1(out.depth, out.pdd.(out.engines{e}), d));
        end
        matRad_cfg.dispInfo('\n');
    end
    matRad_cfg.dispInfo('%8s', 'time [s]');
    for e = 1:numel(out.engines)
        matRad_cfg.dispInfo(' %9.1f', out.time.(out.engines{e}));
    end
    matRad_cfg.dispInfo('\nRoughness (rms 2nd difference / mean, per region):\n');
    for e = 1:numel(out.engines)
        matRad_cfg.dispInfo('%8s %s\n', out.engines{e}, mat2str(out.roughness.(out.engines{e}), 3));
    end
    matRad_printProfileMetrics(out);
    if ~isempty(out.slab)
        matRad_printLungComparison(out);
    end
end
matRad_printFieldSizes(results.fields);
conv = results.convergence;
matRad_cfg.dispInfo('\n=== scatterCutOff convergence (water, 20 x 20 cm) ===\n');
matRad_cfg.dispInfo('%10s %10s %10s %10s\n', 'cut [mm]', 'D(100)', 'D(200)', 'time [s]');
for i = 1:numel(conv.scatterCutOff)
    matRad_cfg.dispInfo('%10g %10.4f %10.4f %10.1f\n', conv.scatterCutOff(i), conv.dose100(i), conv.dose200(i), conv.time(i));
end
end

function matRad_printProfileMetrics(out)
matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispInfo('Lateral profiles (crossings by linear interpolation):\n%8s %8s %10s %12s\n', 'd [mm]', 'engine', 'w50 [mm]', 'p80-20 [mm]');
for d = [15 50 100 200]
    for e = 1:numel(out.engines)
        m = matRad_profileMetrics(out.lateral, out.profiles.(out.engines{e}).(sprintf('d%d', d)));
        matRad_cfg.dispInfo('%8g %8s %10.1f %12.1f\n', d, out.engines{e}, m.width50, m.penumbra);
    end
end
end

function m = matRad_profileMetrics(x, profile)
x = x(:);
profile = profile(:);
c = find(abs(x) == min(abs(x)), 1);
p = profile / mean(profile(max(c - 1, 1):min(c + 1, numel(profile))));
levels = [0.8 0.5 0.2];
right = zeros(1, 3);
left = zeros(1, 3);
for l = 1:3
    right(l) = matRad_profileCrossing(x(c:end), p(c:end), levels(l));
    left(l) = matRad_profileCrossing(flipud(x(1:c)), flipud(p(1:c)), levels(l));
end
m.width50 = right(2) - left(2);
m.penumbra = 0.5 * (abs(right(1) - right(3)) + abs(left(1) - left(3)));
end

function xc = matRad_profileCrossing(x, p, level)
i = find(p < level, 1);
if isempty(i) || i == 1
    xc = NaN;
    return
end
xc = x(i - 1) + (level - p(i - 1)) * (x(i) - x(i - 1)) / (p(i) - p(i - 1));
end

function matRad_printFieldSizes(fs)
matRad_cfg = MatRad_Config.instance();
hasMC = any(strcmp(fs.engines, 'ompMC'));
ref10 = fs.fieldSize == 100;
matRad_cfg.dispInfo('\n=== Field size dependence (water, 100 mm depth, normalized to 10 x 10) ===\n%12s', 'field');
matRad_cfg.dispInfo(' %9s', fs.engines{:});
if hasMC
    matRad_cfg.dispInfo('  CCC/MC SVD/MC  MC 1 sd');
end
matRad_cfg.dispInfo('\n');
for f = 1:numel(fs.fieldSize)
    matRad_cfg.dispInfo('%9g mm', fs.fieldSize(f));
    vals = struct();
    for e = 1:numel(fs.engines)
        name = fs.engines{e};
        vals.(name) = fs.cax.(name)(f) / fs.cax.(name)(ref10);
        matRad_cfg.dispInfo(' %9.4f', vals.(name));
    end
    if hasMC
        matRad_cfg.dispInfo('  %6.3f %6.3f  %5.1f%%', vals.CCC / vals.ompMC, vals.SVDPB / vals.ompMC, 100 * fs.caxRelStd.ompMC(f));
    end
    matRad_cfg.dispInfo('\n');
end
matRad_cfg.dispInfo('%12s', 'out 30 mm');
vals = struct();
for e = 1:numel(fs.engines)
    name = fs.engines{e};
    vals.(name) = fs.out.(name) / fs.cax.(name)(ref10);
    matRad_cfg.dispInfo(' %9.4f', vals.(name));
end
if hasMC
    matRad_cfg.dispInfo('  %6.3f %6.3f  %5.1f%%', vals.CCC / vals.ompMC, vals.SVDPB / vals.ompMC, 100 * fs.outRelStd.ompMC);
end
matRad_cfg.dispInfo('\n%12s', 'time [s]');
for e = 1:numel(fs.engines)
    matRad_cfg.dispInfo(' %9.0f', sum(fs.time.(fs.engines{e})));
end
matRad_cfg.dispInfo('\n');
end

function matRad_printLungComparison(out)
matRad_cfg = MatRad_Config.instance();
points = [out.slab(1) + 25, mean(out.slab(1:2)), out.slab(2) - 25, out.slab(2) - 5, out.slab(2) + 5, out.slab(2) + 25];
labels = {'lung prox', 'lung mid', 'lung dist', 'interf -5', 'interf +5', 'water +25'};
ref = '';
if any(strcmp(out.engines, 'ompMC'))
    ref = 'ompMC';
end
matRad_cfg.dispInfo('Lung and distal interface (%g x %g mm field):\n%10s %8s', out.fieldSize, out.fieldSize, 'point', 'd [mm]');
matRad_cfg.dispInfo(' %9s', out.engines{:});
if ~isempty(ref)
    matRad_cfg.dispInfo('  CCC/MC SVD/MC');
end
matRad_cfg.dispInfo('\n');
for i = 1:numel(points)
    matRad_cfg.dispInfo('%10s %8g', labels{i}, points(i));
    vals = struct();
    for e = 1:numel(out.engines)
        vals.(out.engines{e}) = interp1(out.depth, out.pdd.(out.engines{e}), points(i));
        matRad_cfg.dispInfo(' %9.4f', vals.(out.engines{e}));
    end
    if ~isempty(ref)
        matRad_cfg.dispInfo('  %6.3f %6.3f  (MC 1 sd %.1f%%)', vals.CCC / vals.(ref), vals.SVDPB / vals.(ref), ...
                            100 * interp1(out.depth, out.relStd.(ref), points(i)));
    end
    matRad_cfg.dispInfo('\n');
end
end

function matRad_plotValidation(results)
cases = {'water', 'slab10', 'slab3'};
for c = 1:numel(cases)
    out = results.(cases{c});
    figure('Name', sprintf('CCC validation: %s', cases{c}));
    subplot(1, 2, 1);
    hold on;
    for e = 1:numel(out.engines)
        plot(out.depth, out.pdd.(out.engines{e}), 'DisplayName', out.engines{e});
    end
    xlabel('depth [mm]');
    ylabel('dose');
    title('central axis');
    legend('show');
    subplot(1, 2, 2);
    hold on;
    for e = 1:numel(out.engines)
        prof = out.profiles.(out.engines{e});
        for d = [15 50 100 200]
            plot(out.lateral, prof.(sprintf('d%d', d)), 'DisplayName', sprintf('%s %d mm', out.engines{e}, d));
        end
    end
    xlabel('lateral position [mm]');
    ylabel('dose');
    title('profiles');
end
end
