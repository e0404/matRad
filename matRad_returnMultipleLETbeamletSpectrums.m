function [PhysDose,LET] = matRad_returnMultipleLETbeamletSpectrums(definedEnd,index,ct,ctCube,dij,resultGUI,bins,displayfigures,add)
% matRad histogram plot for multiple LET beamlet Spectrums for different Voxels 
% 
% call
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnMultipleLETbeamletSpectrums(index,ct,ctCube,dij,resultGUI)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnMultipleLETbeamletSpectrums(index,ct,ctCube,dij,resultGUI,displayfigures)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnMultipleLETbeamletSpectrums(index,ct,ctCube,dij,resultGUI,bins)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnMultipleLETbeamletSpectrums(index,ct,ctCube,dij,resultGUI,bins,displayfigures)
%
% input
%   definedEnd:             
%   index:                  voxel coordinates from the cube index in matRadGUI in [y x z]
%   ct:                     matRad ct struct
%   ctCube:                 matRad ct.cube in ct struct
%   dij:                    matRad dij struct
%   resultGUI:              matRad resultGUI struct
%   bins:                   (optional) specifies the number of bins
%   displayfigures:         (optional) displays the figures
%   add:                    add an array [y x z] to go step by step to the next voxel
%
% output
%   LETbeamletSpectrum: histogram plot for LET beamlet spectrum
%   PhysDose_LET:       bivariate histogram plot: dose distribution and associated LET
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LET = [];
PhysDose = [];
    for i = 1:definedEnd                % several spectrums for multiple voxels
        
        [wPhysDose,wLET,~,~] = matRad_plotLETbeamletSpectrumInVoxel(index, ct, ctCube, dij, resultGUI, bins, displayfigures);
        PhysDose = [PhysDose wPhysDose];
        LET = [LET wLET];
        index = index + add;
        
    end
end
