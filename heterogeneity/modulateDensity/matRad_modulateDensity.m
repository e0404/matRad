function [ct] = matRad_modulateDensity(ct,cst,pln,Pmod,mode)
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

% setting overlaps
%lungIdx = lungIdx(~ismember(lungIdx,tumorIdx));

% calculate ct cube from cubeHU if not specified
if ~isfield(ct,'cube')
    ct = matRad_calcWaterEqD(ct,pln);
end

if strcmp(mode, 'binomial')
     
    %pLung = 0.26;
    %pLung = 0.4;
     rhoLung = 1.05;
%    rhoLung = 1;
    
    pLung = ct.cube{1}(lungIdx) / rhoLung;
    if any(pLung > 1)
        lungIdx = lungIdx(pLung <= 1);
        pLung = ct.cube{1}(lungIdx) / rhoLung;
    end
    
    d = Pmod/1000 ./ (1-pLung) / rhoLung; % [1] eq.8: Pmod = d*(1-pLung) * rhoLung
    D = ct.resolution.y;
    
    n = round(D./d);
    
    ct.cube{1}(lungIdx) =  matRad_sampleLungBino(n,pLung,rhoLung,length(lungIdx));
    
elseif strcmp(mode, 'poisson')
    
    DensMod = matRad_loadModDist(Pmod);
    
    DensMod(1,1)=0.001;
    DensMod(:,4)=DensMod(:,2);
    
    %normalize density distribution
%     DensMod(:,4) = DensMod(:,4) / sum(DensMod(:,4));
    
    %Connect each density-probability-couple with a number that will later be transformed the HU-Value:
    %The maximum HU-Value of the HU set by Schneider etal is 2995: So the HU of the modulated density must be at least 2995+1=2996 ; This values is prerocessed with the later used RescaleIntercept and RescaleSlope. See also calculation for Threshold_HU_Value_to_double.
%     Min_HU_for_DensMod_in_double=((2995+1000)/1);
    
    %HU as double in row 2:
%     DensMod(:,2) = (Min_HU_for_DensMod_in_double+1:Min_HU_for_DensMod_in_double+length(DensMod))';
    %Real HU values in row 3:
%     DensMod(:,3) = (Min_HU_for_DensMod_in_double+1:Min_HU_for_DensMod_in_double+length(DensMod))'-1000;
    
    % descrete sampling of the density distribution
    P = [0; cumsum(DensMod(:,4))];
    
    samples = discretize(rand(numel(lungIdx),1),P);
    ct.cube{1}(lungIdx) = samples / max(samples);
end

ct.cubeHU{1}(lungIdx) = 1024*(ct.cube{1}(lungIdx)-1);

% plot histogram of the the lung density distribution
% figure, histogram(ct.cube{1}(lungIdx))
end