clc,clear ,close all

load('/Volumes/WS_exFat/TG119/verification/constRBE/matlab.mat')
load('/Volumes/WS_exFat/TG119/verification/constRBE/resultGUI.mat');
load('/Volumes/WS_exFat/TG119/verification/constRBE/resultGUIMDACC.mat');

description = 'MCN_RBE';

resultGUIMDACC = matRad_calcCubes(VarName2,dij,cst,1);

resultGUI      = matRad_getBeamContributions(resultGUI,cst,stf,dij,'physicalDose');
resultGUIMDACC = matRad_getBeamContributions(resultGUIMDACC,cst,stf,dij,'physicalDose');


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

matRad_calcMultipleDVH(data,cst,pln,Name)


resultGUIvar = resultGUI;

data{1} = resultGUI;
Name{1}= 'constRBE';

data{2} = resultGUIvar;
Name{2}= 'MCN RBE';

matRad_calcMultipleDVH(data,cst,pln,Name)




