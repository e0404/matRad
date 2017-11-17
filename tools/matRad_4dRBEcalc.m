function D = matRad_4dRBEcalc(ct, cst, dij, resultGUI, pln, File)

% calculates all variations of 4D and RBE calculations
%  
%
% input
%   ct:          matRads ct struct
%   cst:         matRads cst struct
%   pln:         matRads pln struct
%   dij:         inclusive alpha beta cubes
%   resultGUI    
%
%   
addpath('D:\Matrad\')
addpath('D:\Matrad\4Ddose')

% output D with all dose cubes
%D.name = {'Dopt', 'Drecalc3D', 'Dopt -Drecalc3D', 'Drecalc4Dconst', 'Drecalc4Dvar', 'Drecalc4Dconst-Drecalc4Dvar','Dopt - Drecalc4Dconst', 'Drecalc3D-Drecalc4Dvar', 'Dopt-Drecalc4Dvar'};
D.name = {'(A)', '(B)', '(A)-(B)', '(C)', '(D)', '(C)-(D)','(A)-(C)', '(B)-(D)', '(A)-(D)'};

D.isolines = {1, 1, 0, 1, 1, 0, 0, 0, 0};
D.fractions = pln.numOfFractions;  

%Pln info
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; 
pln.gantryAngles    = [210]; 
pln.couchAngles     = [0]; 
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';    
pln.numOfFractions         = 30;
pln.runSequencing          = false; 
pln.runDAO                 = false; 
pln.machine                = 'HIT'; 
pln.minNrParticles         = 500000;
pln.LongitudialSpotSpacing = 5;      
pln.calcLET                = false; 

pln.bioOptimization = 'MCN_RBExD';
pln = matRad_getBioModel(pln);
[cst,pln] = matRad_setPlanUncertainties(ct,cst,pln);

%calc dose matrix
stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst,false);

%opt for const RBE
pln.bioOptimization = 'const_RBExD';
pln = matRad_getBioModel(pln);
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst);  
Dopt = resultGUI.RBExDose;
matRad_export_HITXMLPlan_modified('PlanRBE',  pln, stf, resultGUI, 'stfMode')  

%4D for const RBE
[resultGUI, ~] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI,  File);

Drecalc4Dconst = resultGUI.accRBExDose;

%4D for variable RBE
pln.bioOptimization = 'MCN_RBExD';
pln = matRad_getBioModel(pln);
[cst,pln] = matRad_setPlanUncertainties(ct,cst,pln);
[resultGUI, ~] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI,  File);
Drecalc4Dvar = resultGUI.accRBExDose;

%recalc for variable RBE
Drecalc3D  = matRad_calcMcNRBExD(dij, cst, resultGUI);

D.data = {Dopt, Drecalc3D, Dopt-Drecalc3D, Drecalc4Dconst, Drecalc4Dvar, Drecalc4Dconst-Drecalc4Dvar, Dopt - Drecalc4Dconst, Drecalc3D-Drecalc4Dvar, Dopt-Drecalc4Dvar};






