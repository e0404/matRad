%Set up a phantom
builder = matRad_PhantomBuilder([51, 51, 51],[2 2 2],1);

builder.addBoxOAR('PHANTOM', [100 102 100], 'HU', 0, 'objective',DoseObjectives.matRad_SquaredOverdosing(1,0.5));
builder.addBoxOAR('CALIB_POINT', [1 1 1], 'offset', [0 -1 0], 'HU', 0, 'objective',DoseObjectives.matRad_SquaredOverdosing(10,0.5));
builder.addBoxTarget('TARGET_FIELD',[25 1 25], 'offset', [0 25 0], 'HU', 0, 'objective',DoseObjectives.matRad_SquaredDeviation(100,1.0));

[ct, cst] = builder.getctcst();

ct.cube{1} = ones(size(ct.cubeHU{1}));

%% Setup Plan with one beam
pln.radiationMode = 'photons';
pln.machine = 'Generic';
pln.numOfFractions = 1;
pln.propStf.gantryAngles = 0;
pln.propStf.couchAngles = 0;
pln.propStf.addMargin = 0;
pln.propStf.isoCenter = matRad_getIsoCenter(cst,ct) - [0 1 0];
pln.propDoseCalc.doseGrid.resolution = ct.resolution;
pln.propDoseCalc.intConvResolution = 0.5;
pln.propDoseCalc.useCustomPrimaryPhotonFluence = false;
stf = matRad_generateStf(ct,cst,pln);

%% Compute Dose for open field
resultGUI = matRad_calcDoseForward(ct,cst,stf,pln,ones(stf.totalNumOfBixels,1));

%% Compute water equivalent path / profile depth
rt = matRad_RayTracerSiddon(ct.cube, ct);
[alphas, l, rho, d12, ix] = rt.traceRay(stf.isoCenter,stf.sourcePoint,-stf.sourcePoint);

radDepth = cumsum(l .* rho{1}) - l./2;
dose = resultGUI.physicalDose(ix);

calibDose = interp1(radDepth,dose,50);

stf = matRad_computeSSD(stf,ct);
rayPos_bev = vertcat(stf.ray.rayPos_bev);
centralRayIx = find(all(rayPos_bev == [0 0 0],2));

ssd = stf.ray(centralRayIx).SSD;
fprintf('Calibration dose at 50 mm depth in water with SSD = %f cm: %f (should be 1.0).\n', ssd, calibDose);

figure; plot(radDepth,dose);


%%
%figure; 
matRadGUI;


