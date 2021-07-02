function [ct] = matRad_modulateDensity(ct,cst,pln,Pmod,mode,continuous,voxelsize)
% matRad density modulation function
%
% call
%   ct = matRad_modulateDensity(ct,cst,Pmod,mode)
%
% input
%   ct:             ct struct
%   cst:            matRad cst struct
%   Pmod:           Modulation power according to which the modulation will
%                   be created
%   mode:           mode for density modulation ('binominal','poisson')
%                   note: poisson only available for Pmod = 250,450 and 800
%
% output
%   ct:             ct struct with modulated density cube
%
% References
%   [1] https://iopscience.iop.org/article/10.1088/1361-6560/aa641f
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global matRad_cfg;
matRad_cfg =  MatRad_Config.instance();

if nargin < 5
    mode = 'binomial';
    continuous = false;
end

% get all unique tumor indices from PTV segmentations
% idx = strcmp(cst(:,2),'PTV');
% tumorIdx = [cst{idx,4}];
% tumorIdx = unique(vertcat(tumorIdx{:}));

% get all unique lung indices from lung segmentations
idx = contains(cst(:,2),'Lung');
if sum(idx)==0
    matRad_cfg.dispError('No lung segmentation found in cst.\n');
end
lungIdx = [cst{idx,4}];
lungIdx = unique(vertcat(lungIdx{:}));


idxTemp = zeros(length(10:voxelsize:59),length(1:voxelsize:100),length(1:voxelsize:100));
counter = 2;
for e = 1:size(idxTemp,1)
    for r = 1:size(idxTemp,2)
        for t = 1:size(idxTemp,3)
            idxTemp(e,r,t) = counter;
            counter = counter +1;
        end
    end
end

% setting overlaps
%lungIdx = lungIdx(~ismember(lungIdx,tumorIdx));

% calculate ct cube from cubeHU if not specified
if ~isfield(ct,'cube')
    ct = matRad_calcWaterEqD(ct,pln);
end

if strcmp(mode, 'binomial')
    rhoLung = 1.05;
    
    pLung = ct.cube{1}(lungIdx) / rhoLung;
    if any(pLung > 1)
        lungIdx = lungIdx(pLung <= 1);
        pLung = ct.cube{1}(lungIdx) / rhoLung;
    end
      
    d = Pmod/1000 ./ (1-pLung) / rhoLung; % [1] eq.8: Pmod = d*(1-pLung) * rhoLung
    
    if (pln.propStf.gantryAngles == 0 || pln.propStf.couchAngles == 0) && pln.propStf.numOfBeams == 1 && (isfield(pln.propHeterogeneity,'angleCorrection') && pln.propHeterogeneity.angleCorrection)
        angleCorrection = 1/cosd(pln.propStf.gantryAngles);
    else
        angleCorrection = 1;
    end
    D = ct.resolution.y*voxelsize*angleCorrection;
    
    if continuous
       n = D./d;
       tooSmall = n>1;
    else
       n = round(D./d);
       tooSmall = n>1;
    end
    
    % Don't modulate voxel with less than 1 substructures
    lungIdx = lungIdx(tooSmall);
    pLung = pLung(tooSmall);
    n = n(tooSmall);
    
%     discrete = matRad_sampleLungBino(n,pLung,rhoLung,length(lungIdx));
%     cont = matRad_sampleLungBino(n,pLung,rhoLung,length(lungIdx),1);
%     figure
%     map = brewermap(2,'Set1'); 
%     h1 = histogram(cont,'facecolor',map(1,:),'facealpha',.5);
%     figure
%     h2 = histogram(discrete,'facecolor',map(2,:),'facealpha',.5);   
    lungVoxNew = numel(unique(idxTemp(:)));
    samplesTmp = matRad_sampleLungBino(unique(n),unique(pLung),rhoLung,lungVoxNew,continuous);
    samplesTmp = reshape(samplesTmp,[size(idxTemp,1),size(idxTemp,2),size(idxTemp,3)]);
    samples = imresize3(samplesTmp, voxelsize, 'nearest');
    ct.cube{1}(lungIdx) = samples(1:50,1:100,1:100);
    ct.cubeHU{1}(lungIdx) = 1024*(ct.cube{1}(lungIdx)-1);
        
elseif strcmp(mode, 'poisson')
    
    DensMod = matRad_loadModDist(Pmod);
    
    DensMod(1,1)=0.001;
    DensMod(:,4)=DensMod(:,2);
    
    %normalize density distribution
%     DensMod(:,4) = DensMod(:,4) / sum(DensMod(:,4));
    
    %Connect each density-probability-couple with a number that will later be transformed the HU-Value:
    %The maximum HU-Value of the HU set by Schneider etal is 2995: So the HU of the modulated density must be at least 2995+1=2996 ; This values is prerocessed with the later used RescaleIntercept and RescaleSlope. See also calculation for Threshold_HU_Value_to_double.
    Min_HU_for_DensMod_in_double=((2995+1000)/1);
    
    %HU as double in row 2:
    DensMod(:,2) = (Min_HU_for_DensMod_in_double+1:Min_HU_for_DensMod_in_double+length(DensMod))';
    %Real HU values in row 3:
    DensMod(:,3) = (Min_HU_for_DensMod_in_double+1:Min_HU_for_DensMod_in_double+length(DensMod))'-1000;
    

    Z_mod = ct.cubeHU{1};

    numOfSamples = 50*length(lungIdx);
    
    zz1 = ceil(rand(numOfSamples,1)*length(DensMod));
    zz2 = rand(numOfSamples,1);
    ix = zz2 <= DensMod(zz1,4);
    newDist = zz1(ix);
    Z_mod(lungIdx) = DensMod(newDist(1:numel(lungIdx)),3);
    
    ct.cubeHU{1} = Z_mod;
    
%     ct.cube{1}(lungIdx) = (ct.cubeHU{1}(lungIdx)+1)/1024;

%     % descrete sampling of the density distribution
%     P = [0; cumsum(DensMod(:,4))];  
%     samples = discretize(rand(numel(lungIdx),1),P);
%     ct.cube{1}(lungIdx) = samples / max(samples);
end

if strcmp(pln.propHeterogeneity.mode,'TOPAS') && strcmp(mode, 'binomial')
    lung = ct.cube{1}(lungIdx);
    sampledDensities = unique(lung);
    if strcmp(pln.propMC.materialConverter,'HUToWaterSchneider_Lung')
        ct.cubeHU{1}(lungIdx(lung == 1.05)) = 0;
        ct.cubeHU{1}(lungIdx(lung == 0)) = -999;
    elseif strcmp(pln.propMC.materialConverter,'HUToWaterSchneider_mod')
        ct.cubeHU{1}(lungIdx(lung == 1.05)) = 3020;
        ct.cubeHU{1}(lungIdx(lung == 0)) = 2997;
    elseif strcmp(pln.propMC.materialConverter,'HUToWaterSchneider_custom')
        for i = 1:numel(sampledDensities)
            ct.cubeHU{1}(lungIdx(lung == sampledDensities(i))) = 2995 + i;
        end
    end
    
    sampledDensities(1) = 0.001225;
    ct.sampledDensities = sampledDensities;
end

ct.modulated = 1;
% plot histogram of the the lung density distribution
% figure, histogram(ct.cube{1}(lungIdx))
end