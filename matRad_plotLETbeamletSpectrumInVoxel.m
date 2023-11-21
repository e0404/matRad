function [wPhysDose,LET] = matRad_plotLETbeamletSpectrumInVoxel(index, ct, ctCube, dij, resultGUI, bins, displayfigures)
% matRad histogram plot for LET beamlet Spectrum for a certain Voxel &
% bivariate histogram plot to see the dose distribution and associated LET
% 
% call
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_plotLETbeamletSpectrumInVoxel(index,ct,ctCube,dij,resultGUI)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_plotLETbeamletSpectrumInVoxel(index,ct,ctCube,dij,resultGUI,displayfigures)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_plotLETbeamletSpectrumInVoxel(index,ct,ctCube,dij,resultGUI,bins)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_plotLETbeamletSpectrumInVoxel(index,ct,ctCube,dij,resultGUI,bins,displayfigures)
%
% input
%   index:                  voxel coordinates from the cube index in matRadGUI in [y x z]
%   ct:                     matRad ct struct
%   ctCube:                 matRad ct.cube in ct struct
%   dij:                    matRad dij struct
%   resultGUI:              matRad resultGUI struct
%   bins:                   (optional) specifies the number of bins
%   displayfigures:         (optional) displays the figures
%
% output
%   LETbeamletSpectrum: histogram plot for LET beamlet spectrum
%   PhysDose_LET:       bivariate histogram plot: dose distribution and associated LET
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ind = sub2ind(size(ct.cube{ctCube}),index(1),index(2),index(3)); % choose one voxel
    rowOfmLETDose = full(dij.mLETDose{ctCube}(ind,:));               % returns row in the mLETDose matrix for the voxel
    rowDijPhysicalDose = full(dij.physicalDose{ctCube}(ind,:));      % returns row in the physicalDose matrix for the voxel
    
    wLETDose = rowOfmLETDose.*resultGUI.w';                          % calculates weighted dose averaged LET
    wPhysDose = rowDijPhysicalDose.*resultGUI.w';                    % calculates weighted physical dose
     
    LET = wLETDose./wPhysDose;                                       % calculates pure LET

    if ~exist('displayfigures','var') || isempty(displayfigures)     % only if displayfigures is 1 it will return the LET beamlet Spectrum as histogram
        displayfigures = [];
    elseif displayfigures == 1
        figure
        LETbeamletSpectrum = histogram(LET);
        xlabel('LET'); ylabel('number of bixel')

        if ~exist('bins', 'var') || isempty(bins)                    % here the number of bins is defined
            bins = [];
        else
            Nbins = morebins(LETbeamletSpectrum);
            LETbeamletSpectrum.NumBins = bins;
        end

        figure
        PhysDose_LET = histogram2(LET,wPhysDose);                    % bivariate histogram to see the dose distribution and associated LET
        xlabel('LET'); ylabel('physical Dose'); zlabel('number of bixel')
    end
end
