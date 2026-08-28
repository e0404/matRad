%% Example Photon Treatment Plan with VMAT direct aperture optimization
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% In this example we will show
% (i) how to load patient data into matRad
% (ii) how to input necessary parameters in the pln structure
% (iii) how to setup a photon dose calculation
% (iv) how to inversely optimize fluence directly from command window in MatLab.
% (v) how to apply a sequencing algorithm
% (vi) how to run a VMAT direct aperture optimization
% (vii) how to visually and quantitatively evaluate the result
%
% VMAT is direct aperture optimization on a rotating arc. For the static
% (step-and-shoot) counterpart - including a comparison of the siochi, xia
% and engel leaf sequencers - see matRad_example3_photonsDAO.

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the TG119 patient
% into your workspace
matRad_rc;

load TG119.mat;

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This
% structure requires input from the treatment planner and defines
% the most important cornerstones of your treatment plan.

% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'photons';
pln.machine         = 'Generic';

pln.bioModel = 'none';      % biological RBE model, not interesting for photons
pln.multScen = 'nomScen';   % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'

% beam geometry settings
pln.propStf.generator                = 'PhotonVMAT';
pln.propStf.bixelWidth               = 5;           % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles             = [-180, 180]; % gantry arc anchor points
pln.propStf.couchAngles              = [0, 0];      % couch angle for arcs
% pln.propStf.arcIndex                 = [1 1];     % assign anchor points to arcs (if more than one arc is defined)
pln.propStf.maxGantryAngleSpacing    = 8;           % [deg] / max gantry angle spacing for dose calculation
pln.propStf.maxDAOGantryAngleSpacing = 8;           % [deg] / max gantry angle spacing for DAO
pln.propStf.maxFMOGantryAngleSpacing = 24;          % [deg] / max gantry angle spacing for FMO
pln.propStf.minAperturesPerFMOBeam   = 3;           % Minimum number of apertures sequenced for DAO per fluence-optimized beam.

% The FMO spacing is the parameter that decides how faithfully the arc can be
% sequenced, so it should stay small. Each fluence map is optimized at one
% gantry angle but delivered as minAperturesPerFMOBeam separate apertures
% spread over the FMO sector - the sum reproduces the optimized fluence, yet
% every individual control point receives only one fragment of it. The wider
% that sector, the more the delivered dose clusters into a few entrance
% directions. Measured on this case as the peak-to-mean of the out-of-target
% dose binned by gantry angle, the excess introduced by sequencing is +2% at
% 24 deg, +30% at 45 deg and +38% at 72 deg - and direct aperture
% optimization amplifies whatever it inherits rather than repairing it.
% Target coverage and maximum dose look healthy in every case, so this shows
% up only in the spatial dose distribution.

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
pln.propDoseCalc.precision = 'single';

% sequencing settings
pln.propSeq.sequencer          = 'siochi'; % the only sequencer with VMAT support

% numLevels sets the number of stratification levels for decomposing each
% FMO fluence map into apertures. For VMAT it is only a starting value: the
% sequencer re-decomposes at increasing level counts until every FMO beam
% yields at least as many apertures as it has DAO child angles, since each
% DAO control point receives exactly one aperture.
pln.propSeq.numLevels          = 10;

pln.propSeq.continuousAperture = false;  % interpolate leaf positions between DAO control points (dynamic delivery)
pln.propSeq.preconditioner     = true;   % apply Jacobi preconditioning to the aperture weights

% Smoothing of each FMO fluence map before it is decomposed into apertures
% ('gaussian' or 'none'). The blurred field rim is what lets the sequencer
% produce one aperture per DAO control point, at the price of the sequenced
% fluence no longer reproducing the optimized one - which is why it is
% applied for VMAT only. 'none' is only usable for fluences that are already
% modulated enough to stratify into as many apertures as there are DAO
% control points per FMO beam.
pln.propSeq.arcFluenceSmoothing = 'gaussian';

% The stratification generates more apertures per FMO beam than the arc can
% deliver, so all but numOfChildren of them are discarded. 'doseAreaProduct'
% keeps the ones with the largest open area and rescales them by a common
% factor; 'leastSquares' instead keeps the apertures that best reconstruct
% the optimized fluence and refits their weights to it, preserving the total
% dose-area product. On this example the latter starts the DAO closer to the
% fluence optimization result, mainly by producing a less hot plan.
pln.propSeq.apertureSelection   = 'leastSquares';

% optimization settings
pln.propOpt.quantityOpt         = 'physicalDose';   % Quantity to optimizer (could also be RBExDose, BED, effect)
pln.propOpt.optimizer           = 'IPOPT';          % We can also utilize 'fmincon' from Matlab's optimization toolbox
pln.propOpt.runVMAT             = true;             % put fluence optimization, sequencing and DAO into VMAT (arc) mode

% add a fluence objective for smoothing for smoother fluences
pln.propOpt.fluenceObjectives   = {FluenceObjectives.matRad_FluenceVariance(5)};

% optionally, the final plan can be scaled such that the target D95 matches
% the prescription:
% pln.propOpt.scaleToPrescription  = true;
% pln.propOpt.prescribedDose       = 60; % [Gy] over all fractions
% pln.propOpt.prescriptionStructIx = find(strcmp(cst(:, 3), 'TARGET'));

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct, cst, pln);

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence
% matrices for unit beamlet intensities. Having dose influences available
% allows for subsequent inverse optimization.
dij = matRad_calcDoseInfluence(ct, cst, stf, pln);

%% Inverse Planning for IMRT
% The goal of the fluence optimization is to find a set of beamlet weights
% which yield the best possible dose distribution according to the
% predefined clinical objectives and constraints underlying the radiation
% treatment. In VMAT, FMO is done only at the FMO angle subset of the arc;
% the dose influence matrix carries this bookkeeping (dij.isFMOBeam) from
% the steering information. Once the optimization has finished, trigger
% once the GUI to visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij, cst, pln);
matRadGUI;

%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in
% order to modulate the intensity of the beams with multiple static
% segments, so that translates each intensity map into a set of deliverable
% aperture shapes. The fluence map at each angle in the initGantryAngles
% set is sequenced, with the resulting apertures spread to neighbouring
% angles from the optGantryAngles set.
resultGUI = matRad_sequencing(resultGUI, stf, pln, dij);

%% DAO - Direct Aperture Optimization
% The Direct Aperture Optimization is an optimization approach where we
% directly optimize aperture shapes and weights at the angles in the
% optGantryAngles set.  The gantry angle speed, leaf speed, and MU rate are
% constrained by the min and max values specified by the user.
resultGUI = matRad_directApertureOptimization(dij, cst, resultGUI.apertureInfo, resultGUI, pln);

%% Aperture visualization
% Use a matRad function to visualize the result.

visLeafCoordinate = 'physical'; % 'leafNum' to show leaf numbers
visMode = 'animate'; % other options 'perBeam', 'trajectory', 'interactive'
matRad_visApertureInfo(resultGUI.apertureInfo, visMode, 'leafCoordinate', visLeafCoordinate);

%% Indicator Calculation and display of DVH and QI
resultGUI = matRad_planAnalysis(resultGUI, ct, cst, stf, pln);

%% Calculate delivery metrics
% Total MU and delivery time of the plan plus per-control-point MU rate,
% gantry speed and leaf speeds are collected in resultGUI.deliveryMetrics.
resultGUI = matRad_calcDeliveryMetrics(resultGUI, stf);

%% Dose recalculation on a finer arc
% The plan was optimized at the (coarse) DAO control points. To check the
% dose that is actually delivered along the arc, resample the optimized
% apertures onto a finer gantry-angle grid and forward-calculate the dose.
refinedSpacing = 2; % [deg] / gantry angle spacing for dose calculation
[stfFine, apertureInfoFine] = matRad_refineApertureArc(ct, cst, pln, resultGUI.apertureInfo, refinedSpacing);
resultGUIfine = matRad_calcDoseForward(ct, cst, stfFine, pln, apertureInfoFine.bixelWeights);
resultGUI.physicalDose_fine = resultGUIfine.physicalDose;

% Instead of recomputing the dose, an already available dij can also be
% recycled by redirecting the fine beams onto the original DAO angles: see
% the 'reuseDij' option of matRad_refineApertureArc.
