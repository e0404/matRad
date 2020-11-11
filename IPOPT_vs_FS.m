load TG119.mat

pln.radiationMode = 'photons';  
pln.machine       = 'Generic';
pln.propOpt.bioOptimization = 'none';    
pln.numOfFractions         = 1;
pln.propStf.gantryAngles   = [0:72:359];
%pln.propStf.gantryAngles   = [0];
pln.propStf.couchAngles    = zeros(1, numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]
pln.propOpt.runSequencing = 1;
pln.propOpt.runDAO        = 0;
stf                      = matRad_generateStf(ct,cst,pln);
dij = matRad_calcPhotonDose(ct,stf,pln,cst);
pln.propOpt.tol_obj = 1e-20;
pln.propOpt.tol_violation = 1e-20;
pln.propOpt.accepted_violation = 1e-20;

cst{1, 5}.visibleColor = [0 0.4470 0.7410];
cst{2, 5}.visibleColor = [0.9290 0.6940 0.1250];
cst{3, 5}.visibleColor = [0.3010 0.7450 0.9330];

cst = cst(:, 1:5);
cst{1, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(100, 15)
cst{2, 6}{1}=DoseObjectives.matRad_SquaredDeviation(1000, 50);
cst{3, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(30, 25);

cst{2, 6}{1}=DoseConstraints.matRad_MinMaxDose(49, 51, 1, 1000);
%cst{1, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 20, 1, 100);
%cst{3, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 30, 1, 30);

%% IPOPT
pln.propOpt.optimizer = "IPOPT";
%pln.propOpt.feasibility_seeker = "AMS_sim";
pln.propOpt.max_iter = 500;
pln.propOpt.max_time = 3600;
%pln.propOpt.lambda = 1;
%pln.propOpt.weighted = true;
tic;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
time = toc;


%% Super AMS_sim
pln.propOpt.optimizer = "Superization";
pln.propOpt.feasibility_seeker = "AMS_sim";
pln.propOpt.max_iter = 500;
pln.propOpt.max_time = 3600;
pln.propOpt.lambda = 1;
pln.propOpt.weighted = true;

tic;
resultGUI_super = matRad_fluenceOptimization(dij,cst,pln);
time = toc;

%% 

plane      = 3;
doseWindow = [-30 30];
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

absDiffCube = resultGUI.physicalDose-resultGUI_super.physicalDose;
figure,title( 'IPOPT plan - Superiorization plan')
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
%%

plane      = 3;
doseWindow = [0 max([resultGUI.physicalDose(:); resultGUI_super.physicalDose(:)])];

figure,title('IPOPT plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('Superiorization plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_super.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);

%% Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
[dvh,qi]               = matRad_indicatorWrapper(cst,pln,resultGUI);
[dvh_super,qi_super] = matRad_indicatorWrapper(cst,pln,resultGUI_super);

%% Show DVHs
figure()
matRad_showDVH(dvh,cst,pln,1)
hold on
matRad_showDVH(dvh_super,cst,pln,2)

