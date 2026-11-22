%% Example: Photon treatment plan with collapsed cone convolution dose calculation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show
% (i) how to compile the collapsed cone transport mex file (optional)
% (ii) how to set up a phantom with a lung insert using the phantom builder
% (iii) how to calculate the dose influence matrix of a photon plan with the
% collapsed cone convolution engine
% (iv) how to recalculate the optimized plan with the SVD pencil beam
% engine and compare both dose distributions

%% Set matRad runtime configuration
% Let's start with a clean environment and the matRad configuration object,
% which we use for the console output of this example.
matRad_rc;
matRad_cfg = MatRad_Config.instance();

%% Compile the transport (optional)
% The collapsed cone engine transports the energy released by each pencil
% beam (TERMA) along a set of cone directions. This transport is implemented
% as C++ mex file. If the mex file is not compiled, the engine automatically
% uses a MATLAB implementation of the same algorithm, which gives the same
% result but is slower. The mex file only needs to be compiled once.
%
% With MATLAB on Linux and a compiler newer than the one supported by your
% MATLAB release, use matRad_compileCollapsedConeMex([], true, true) to link
% the C++ runtime statically.
if exist('matRad_collapsedConeMex', 'file') ~= 3
    try
        matRad_compileCollapsedConeMex();
    catch err
        matRad_cfg.dispWarning('Collapsed cone mex file not compiled: %s', err.message);
    end
end

%% Phantom with a lung insert
% We create a water phantom of 80 x 80 x 50 voxels of 3 mm with the
% phantom builder. Inside, a box of lung tissue (-750 HU, i.e. a relative
% electron density of about 0.25) contains a spherical tumor of 6 voxels
% radius. Each volume gets an objective for the optimization: a prescribed
% dose of 60 Gy for the tumor and overdosing penalties for the lung and the
% remaining body.
%
% The order matters: the volumes defined first have priority where they
% overlap, so the tumor is defined before the lung and the lung before the
% body.
builder = matRad_PhantomBuilder([80 80 50], [3 3 3], 1);
objTumor = struct(DoseObjectives.matRad_SquaredDeviation(800, 60));
objLung = struct(DoseObjectives.matRad_SquaredOverdosing(100, 20));
objBody = struct(DoseObjectives.matRad_SquaredOverdosing(10, 30));
builder.addSphericalTarget('Tumor', 6, 'objectives', objTumor, 'HU', 0);
builder.addBoxOAR('Lung', [36 40 30], 'objectives', objLung, 'HU', -750);
builder.addBoxOAR('Body', [76 76 46], 'objectives', objBody, 'HU', 0);
[ct, cst] = builder.getctcst();

%% Treatment plan
% We define a photon plan with five equidistant beams and 10 mm bixels,
% using the generic 6 MV machine.
pln.radiationMode = 'photons';
pln.machine = 'Generic';
pln.numOfFractions = 30;
pln.bioModel = 'none';
pln.multScen = 'nomScen';
pln.propStf.gantryAngles = 0:72:359;
pln.propStf.couchAngles = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth = 10;
pln.propStf.isoCenter = matRad_getIsoCenter(cst, ct, 0);
pln.propOpt.runDAO = false;
pln.propSeq.runSequencing = false;

%%
% The collapsed cone engine is selected with its short name 'CCC'. The
% number of polar and azimuthal cones and the scatter cut-off, i.e. the
% lateral margin around each pencil beam in which the scattered energy is
% transported, can be set in pln.propDoseCalc (the values below are the
% defaults). More cones and a larger cut-off increase the accuracy and the
% computation time. All options are described in the pln documentation.
pln.propDoseCalc.engine = 'CCC';
pln.propDoseCalc.numPolarAngles = 12;
pln.propDoseCalc.numAzimuthalAngles = 8;
pln.propDoseCalc.scatterCutOff = 150;
pln.propDoseCalc.doseGrid.resolution.x = 3;
pln.propDoseCalc.doseGrid.resolution.y = 3;
pln.propDoseCalc.doseGrid.resolution.z = 3;

%% Generate beam geometry
% The steering information (beams, rays and bixels) is generated as for any
% other photon plan.
stf = matRad_generateStf(ct, cst, pln);

%% Dose influence matrix with the collapsed cone engine
% The dose influence matrix is computed bixel by bixel. For each bixel the
% engine computes the TERMA along the ray, transports it with the collapsed
% cones through the relative electron density of the CT, and stores the
% dose in the dij struct exactly like the pencil beam engines. The dose is
% normalized to the depth dose of the machine data in water.
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

%% Fluence optimization
% The optimization does not depend on the dose engine: it only uses the
% dose influence matrix.
resultGUI = matRad_fluenceOptimization(dij, cst, pln);

%% Recalculation with the SVD pencil beam engine
% To see the effect of the dose engine, we recalculate the optimized plan
% with the singular value decomposed pencil beam engine, using the same
% bixel weights and dose grid.
plnSVD = pln;
plnSVD.propDoseCalc = struct('engine', 'SVDPB', 'doseGrid', pln.propDoseCalc.doseGrid);
resultSVD = matRad_calcDoseForward(ct, cst, stf, plnSVD, resultGUI.w);

%% Comparison
% The pencil beam algorithm scales the dose with the radiological depth
% only. For small fields in and behind low density tissue the lateral
% electron equilibrium is not reached, which the pencil beam algorithm does
% not model: it overestimates the dose in the lung and in a tumor
% surrounded by lung. The collapsed cone algorithm transports the scattered
% energy through the actual densities and shows the broader penumbra and
% the lower dose in the lung. matRad_compareDose shows both dose
% distributions, their difference and the dose volume histograms.
matRad_compareDose(resultGUI.physicalDose, resultSVD.physicalDose, ct, cst, [1 0 1], 'on', pln);

%%
% Finally we print the mean dose of each structure for both engines.
for i = 1:size(cst, 1)
    voxels = cst{i, 4}{1};
    matRad_cfg.dispInfo('%-6s mean dose CCC %6.2f Gy, SVDPB %6.2f Gy\n', cst{i, 2}, ...
                        mean(resultGUI.physicalDose(voxels)), mean(resultSVD.physicalDose(voxels)));
end
