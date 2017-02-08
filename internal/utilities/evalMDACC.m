clc,clear ,close all

load('/Volumes/WS_exFat/TG119/verification/constRBE/matlab.mat')
load('/Volumes/WS_exFat/TG119/verification/constRBE/resultGUI.mat');
load('/Volumes/WS_exFat/TG119/verification/constRBE/resultGUIMDACC.mat');


matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.RBExDose,resultGUIMDACC.RBExDose,{'matRad','MDACC'},65,3)

matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose,resultGUIMDACC.physicalDose,{'matRad','MDACC'},65,3,'varRBE_physDose')
matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.RBExDose,resultGUIMDACC.RBExDose,{'matRad','MDACC'},65,3,'varRBE_RBExD')


matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose_beam_1,resultGUIMDACC.physicalDose_beam_1,{'matRad','MDACC'},65,3,'varRBE_beam1')
matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose_beam_2,resultGUIMDACC.physicalDose_beam_2,{'matRad','MDACC'},65,3,'varRBE_beam2')
matRad_plotTwoDoseCubes(ct,cst,pln,resultGUI.physicalDose_beam_3,resultGUIMDACC.physicalDose_beam_3,{'matRad','MDACC'},65,3,'varRBE_beam3')





data{1} = resultGUI;
Name{1}= 'matRad';

data{2} = resultGUIMDACC;
Name{2}= 'MDACC';

matRad_calcMultipleDVH(data,cst,pln,Name)




