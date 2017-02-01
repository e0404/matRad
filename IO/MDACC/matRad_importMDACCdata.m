function [dij] = matRad_importMDACCdata(patientDir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

extension = 'ASCII';

%ijFileNames = matRad_convertMDACCInfluenceDataToMAT(patientDir,extension);


ijFileNames{1,1} = '//Volumes/WS_exFat/patient1/Dij_1.mat';
ijFileNames{2,1} = '/Volumes/WS_exFat/patient1/Dij_2.mat';
%ijFileNames{3,1} = '/Users/Hans-PeterWieser/Documents/Heidelberg/MDAnderson2016/patient1/LETij_1.mat';
%ijFileNames{4,1} = '/Users/Hans-PeterWieser/Documents/Heidelberg/MDAnderson2016/patient1/LETij_2.mat';

ijFileNames{1,2} = 'Dij';
ijFileNames{2,2} = 'Dij';
ijFileNames{3,2} = 'LETij';
ijFileNames{4,2} = 'LETij';

dij              = matRad_convertMATtoSPARSE(ijFileNames);


end

