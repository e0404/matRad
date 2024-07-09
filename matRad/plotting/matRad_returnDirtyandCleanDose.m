function [highLETphysDose,lowLETphysDose,totalphysDose] = matRad_returnDirtyandCleanDose(index,ct,ctCube,dij,resultGUI,LET_thres,displayfigures,maxDirtyDose,bins)
% matRad function to return histogram plots for the dirty and clean dose share of physical dose
%
% call
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnDirtyandCleanDose(index,ct,ctCube,dij,resultGUI,LET_thres)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnDirtyandCleanDose(index,ct,ctCube,dij,resultGUI,LET_thres,displayfigures)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnDirtyandCleanDose(index,ct,ctCube,dij,resultGUI,LET_thres,displayfigures,maxDirtyDose)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_returnDirtyandCleanDose(index,ct,ctCube,dij,resultGUI,LET_thres,displayfigures,maxDirtyDose,bins)
%
% input
%   index:              voxel coordinates from the cube index in matRadGUI in [y x z]
%   ct:                 matRad ct struct
%   ctCube:             matRad ct.cube in ct struct
%   dij:                matRad dij struct
%   resultGUI:          matRad resultGUI struct
%   LET_thres:          LET threshold: above this = dirty dose, below this = clean dose
%   displayfigures:     (optional) displays the figures
%   maxDirtyDose:       (optional) maximum allowed dirty dose, dirty dose in a voxel exceeding this value gets penalized
%   bins:               (optional) specifies the number of bins
%
% output
%   shareDose: bar plot for comparing dirty, clean and the total dose; each summed up to one value
%   dirtyDose: histogram plot for the dirty dose share of physical dose
%   cleanDose: histogram plot for the clean dose share of physical dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [wPhysDose,LET] = matRad_plotLETbeamletSpectrumInVoxel(index, ct, ctCube, dij, resultGUI,[],[]);
    
    highLET = LET > LET_thres;                                  % defining high LET above a certain threshold
    
    % share of high LET contribution in physical dose
    highLETphysDose = sum(highLET.*wPhysDose);                  % sum of all the physical dose with LET value higher than the threshold
    lowLETphysDose = sum(~highLET.*wPhysDose);                  % sum of all the physical dose with LET value lower than the threshold
    totalphysDose = sum(wPhysDose);                             % sum of the total dose deposit in this voxel
    
    if ~exist('displayfigures','var') || isempty(displayfigures)    % only if displayfigures is 1, it will return a bar plot
        displayfigures = [];
    elseif displayfigures == 1
        figure
        dose = [highLETphysDose lowLETphysDose totalphysDose];
        C2 = {'dirty share of physical dose';'clean share of physical dose';'total physical dose in this voxel'};
        shareDose = bar(dose);                                  % summary of the composition of physical dose
        xticklabels(C2)
    end
    
    dirtydose = wPhysDose(highLET);                             % dirty dose defined as the physical dose deposited by the bixels with high LET
    cleandose = wPhysDose(~highLET);                            % clean dose defined as the physical dose deposited by the bixels with low LET
    
    if ~exist('maxDirtyDose', 'var') || isempty(maxDirtyDose)   % if there is a maximum dose up to which the dirty dose should be considered, remove all the dose values above
        maxDirtyDose = [];
    else
        dirtydose(dirtydose>maxDirtyDose) = [];
    end
    
    if ~exist('displayfigures','var') || isempty(displayfigures)
        displayfigures = [];
    elseif displayfigures == 1
        figure
        dirtyDose = histogram(dirtydose);                       % plot the dirty dose as histogram
    
        if ~exist('bins', 'var') || isempty(bins)
            bins = [];
        else
            Nbins = morebins(dirtyDose);                        % define the number of bins
            dirtyDose.NumBins = bins;
        end
        xlabel('dirty dose'); ylabel('number of bixel')
    
        figure
        cleanDose = histogram(cleandose);                       % plot the clean dose as histogram
    
        if ~exist('bins', 'var') || isempty(bins)
            bins = [];
        else
            Nbins = morebins(cleanDose);                        % define the number of bins
            cleanDose.NumBins = bins;
        end
        xlabel('clean dose'); ylabel('number of bixel')
    end
    
end