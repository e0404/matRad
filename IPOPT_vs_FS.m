load TG119.mat

pln.radiationMode = 'photons';  
pln.machine       = 'Generic';
pln.propOpt.bioOptimization = 'none';    
pln.numOfFractions         = 1;
pln.propStf.gantryAngles   = [0:50:359];
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

%% IPOPT
cst = cst(:, 1:5);
cst{1, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(100, 20)
cst{2, 6}{1}=DoseObjectives.matRad_SquaredDeviation(1000, 50);
cst{3, 6}{1}=DoseObjectives.matRad_SquaredOverdosing(30, 30);

pln.propOpt.optimizer = "IPOPT";
%pln.propOpt.feasibility_seeker = "AMS_sim";
pln.propOpt.max_iter = 100;
pln.propOpt.max_time = 3000;
%pln.propOpt.lambda = 1;
%pln.propOpt.weighted = true;
tic;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
time = toc;


%% Super AMS_sim

cst = cst(:, 1:5);
cst{1, 6}{2}=DoseObjectives.matRad_SquaredOverdosing(100, 20)
cst{2, 6}{2}=DoseObjectives.matRad_SquaredDeviation(1000, 50);
cst{2, 6}{1}=DoseConstraints.matRad_MinMaxDose(49, 51, 1, 1.9);
cst{1, 6}{1}=DoseConstraints.matRad_MinMaxDose(0, 20, 1, 0.9);
pln.propOpt.optimizer = "Superization";
pln.propOpt.feasibility_seeker = "AMS_sim";
pln.propOpt.max_iter = 100;
pln.propOpt.max_time = 3000;
pln.propOpt.lambda = 1;
pln.propOpt.weighted = true;

tic;
resultGUI_coarse = matRad_fluenceOptimization(dij,cst,pln);
time = toc;

%% 

plane      = 3;
doseWindow = [-30 30];
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);

absDiffCube = resultGUI.physicalDose-resultGUI_coarse.physicalDose;
figure,title( 'fine beam spacing plan - coarse beam spacing plan')
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
%%

plane      = 3;
doseWindow = [0 max([resultGUI.physicalDose(:); resultGUI_coarse.physicalDose(:)])];

figure,title('original plan - fine beam spacing')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('modified plan - coarse beam spacing')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_coarse.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);

%% Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
[dvh,qi]               = matRad_indicatorWrapper(cst,pln,resultGUI);
[dvh_coarse,qi_coarse] = matRad_indicatorWrapper(cst,pln,resultGUI_coarse);

figure()
matRad_showDVH(dvh,cst,pln,1)
hold on
matRad_showDVH(dvh_coarse,cst,pln,2)

