function [ resultGUI ] = matRad_getBeamContributions(resultGUI,cst,stf,dij,quantity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tmpResultGUI.(quantity) = resultGUI.(quantity);
tmpResultGUI.physicalDose = resultGUI.physicalDose;

ix = 1;
for i = 1:size(stf,2)
   
   tmpResultGUI.w = zeros(size(resultGUI.w));
   
   tmpResultGUI.w(ix:ix+stf(i).totalNumOfBixels-1) = resultGUI.w(ix:ix+stf(i).totalNumOfBixels-1);
   
   
    BeamContrib = matRad_calcCubes(tmpResultGUI.w,dij,cst,1);
   
    resultGUI.([quantity '_beam_' num2str(i,'%d')]) = BeamContrib.(quantity);
   ix = ix + stf(i).totalNumOfBixels;
end




end

