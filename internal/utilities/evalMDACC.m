clc,clear ,close all

load('/Volumes/WS_exFat/TG119/verification/MCN_RBExD/matlab.mat')
load('/Volumes/WS_exFat/TG119/verification/MCN_RBExD/resultGUI.mat');
load('/Volumes/WS_exFat/TG119/verification/MCN_RBExD/resultGUIMDACC.mat');

description = 'MCN_RBExD_TG119';
quantity    = 'RBExDose';%'physicalDose';

%resultGUIMDACC = matRad_calcCubes(VarName2,dij,cst,1);

resultGUI      = matRad_getBeamContributions(resultGUI,cst,stf,dij,quantity);
resultGUIMDACC = matRad_getBeamContributions(resultGUIMDACC,cst,stf,dij,quantity);


matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.RBExDose,resultGUIMDACC.RBExDose,{'matRad','MDACC'},65,3)

matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.RBExDose,resultGUIMDACC.RBExDose,{'matRad','MDACC'},65,3,description)

matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose,resultGUIMDACC.physicalDose,{'matRad','MDACC'},65,3,description)



matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.RBExDose_beam_1,resultGUIMDACC.RBExDose_beam_1,{'matRad','MDACC'},65,3,[description 'beam1'])
matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.RBExDose_beam_2,resultGUIMDACC.RBExDose_beam_2,{'matRad','MDACC'},65,3,[description 'beam2'])
matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.RBExDose_beam_3,resultGUIMDACC.RBExDose_beam_3,{'matRad','MDACC'},65,3,[description 'beam3'])


matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose_beam_1,resultGUIMDACC.physicalDose_beam_1,{'matRad','MDACC'},65,3,[description 'beam1'])
matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose_beam_2,resultGUIMDACC.physicalDose_beam_2,{'matRad','MDACC'},65,3,[description 'beam2'])
matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose_beam_3,resultGUIMDACC.physicalDose_beam_3,{'matRad','MDACC'},65,3,[description 'beam3'])


data{1} = resultGUI;
Name{1}= 'matRad';

data{2} = resultGUIMDACC;
Name{2}= 'MDACC';

matRad_calcMultipleDVH(data,cst,pln,Name,description,true)
%%

clc,clear ,close all
description = 'compBioModelTG119';
load('/Volumes/WS_exFat/TG119/nominal/constRBE/matlab.mat');

load('/Volumes/WS_exFat/TG119/nominal/constRBE/resultGUI.mat')
constRBE = resultGUI;
load('/Volumes/WS_exFat/TG119/nominal/LSM_RBE/resultGUI.mat')
LSMRBE = resultGUI;
load('/Volumes/WS_exFat/TG119/nominal/MCN_RBE/resultGUI.mat')
MCNRBE = resultGUI;



data{1} = constRBE;
Name{1}= 'constRBE';

data{2} = LSMRBE;
Name{2}= 'LSM';

data{3} = MCNRBE;
Name{3}= 'MCN';


matRad_calcMultipleDVH(data,cst,pln,'RBExDose',Name,description,false)


description = 'compBioModelTG119phys';
matRad_calcMultipleDVH(data,cst,pln,'physicalDose',Name,description,false)



[ constRBE2 ] = matRad_reCalcMCNRBExD(cst,dij,constRBE);
load('/Volumes/WS_exFat/TG119/nominal/LSM_RBE/matlab.mat')
[ LSMRBE2 ]   = matRad_reCalcMCNRBExD(cst,dij,LSMRBE);


data2{1} = constRBE2;
data2{2} = LSMRBE2;
data2{3} = MCNRBE;
description = 'compBioModelTG119recalc';
matRad_calcMultipleDVH(data2,cst,pln,'RBExDose',Name,description,true)





