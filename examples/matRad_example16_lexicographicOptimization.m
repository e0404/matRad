matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
load('TG119.mat');


pln.radiationMode = 'photons';  
pln.machine       = 'Generic';

quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  

%
pln.numOfFractions         = 30;
pln.propStf.gantryAngles   = [0:40:359];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 5;

pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

pln.propSeq.runSequencing = 1;
pln.propOpt.runDAO        = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen');
%%
stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcPhotonDose(ct,stf,pln,cst);

%%
PriorityList1 = matRad_PriorityList1();
PriorityList1.addObjective(1,DoseObjectives.matRad_SquaredDeviation(100,50),4,2);
PriorityList1.addObjective(2,DoseObjectives.matRad_EUD(100,0),10,1);
PriorityList1.addObjective(3,DoseObjectives.matRad_MeanDose(100,0),5,3);

%constraints
%PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,40),1);
%PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(45,57),2);
%PriorityList1.addConstraint(DoseConstraints.matRad_MinMaxDose(0,45),3);
%%
PriorityList1.printPriorityList(cst)
%%
[resultGUIs1,resultGUIs2,cst1,cst2,PriorityList2] =  matRad_2pecOptimization(PriorityList1,dij,cst,pln); 
%%
PriorityList1.printPriorityList(cst)