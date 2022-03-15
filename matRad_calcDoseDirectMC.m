function resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,w,nHistories)
% matRad function to bypass dij calculation for MC dose calculation 
% matRad dose calculation wrapper for MC dose calculation algorithms
% bypassing dij calculation for MC dose calculation algorithms.
% 
% call
%   resultGUI = matRad_calcDoseDirecMC(ct,stf,pln,cst)
%   resultGUI = matRad_calcDoseDirecMC(ct,stf,pln,cst,w)
%   resultGUI = matRad_calcDoseDirectMC(ct,stf,pln,cst,w,nHistories)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   w:          (optional, if no weights available in stf): bixel weight
%               vector
%   nHistories: (optional) number of histories
%
% output
%   resultGUI:  matRad result struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matRad_cfg =  MatRad_Config.instance();

calcDoseDirect = true;

if nargin < 6 || ~exist('nHistories')
  nHistories = matRad_cfg.propMC.direct_defaultHistories;
  matRad_cfg.dispInfo('Using default number of Histories: %d\n',nHistories);
end

% check if weight vector is available, either in function call or in stf - otherwise dose calculation not possible
if ~exist('w','var') && ~isfield([stf.ray],'weight')
     matRad_cfg.dispError('No weight vector available. Please provide w or add info to stf');
end

% copy bixel weight vector into stf struct
if exist('w','var')
    if sum([stf.totalNumOfBixels]) ~= numel(w)
        matRad_cfg.dispError('weighting does not match steering information');
    end
    counter = 0;
    for i = 1:size(stf,2)
        for j = 1:stf(i).numOfRays
            for k = 1:stf(i).numOfBixelsPerRay(j)
                counter = counter + 1;
                stf(i).ray(j).weight(k) = w(counter);
            end
        end
    end
else % weights need to be in stf!
    w = NaN*ones(sum([stf.totalNumOfBixels]),1);
    counter = 0;
    for i = 1:size(stf,2)
        for j = 1:stf(i).numOfRays
            for k = 1:stf(i).numOfBixelsPerRay(j)
                counter = counter + 1;
                w(counter) = stf(i).ray(j).weight(k);
            end
        end
    end    
end       

% dose calculation
if strcmp(pln.radiationMode,'protons')
    engines = {'TOPAS','MCsquare'};
    if ~isfield(pln,'propMC') || ~isfield(pln.propMC,'proton_engine') || ~any(strcmp(pln.propMC.proton_engine,engines))
        matRad_cfg.dispInfo('Using default proton MC engine "%s"\n',matRad_cfg.propMC.default_proton_engine);
        pln.propMC.proton_engine = matRad_cfg.propMC.default_proton_engine;
    end
    
    switch pln.propMC.proton_engine
        case 'MCsquare'
        dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,nHistories,calcDoseDirect);
        case 'TOPAS'
        dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,cst,nHistories,calcDoseDirect);
    end
else
    matRad_cfg.dispError('Forward MC only implemented for protons.');
end


dij.numOfBeams = size(dij.physicalDose{1},2);
dij.beamNum = 1:size(dij.physicalDose{1},2);

% calc resulting dose
if pln.multScen.totNumScen == 1
    % calculate cubes; use uniform weights here, weighting with actual fluence 
    % already performed in dij construction
    resultGUI    = matRad_calcCubes(ones(size(dij.physicalDose{1},2),1),dij,1);
    
% calc individual scenarios    
else    
   Cnt          = 1;
   ixForOpt     = find(~cellfun(@isempty, dij.physicalDose))';
   for i = ixForOpt
      tmpResultGUI = matRad_calcCubes(ones(size(dij.physicalDose{i},2),1),dij,i);
      if i == 1
         resultGUI.([pln.bioParam.quantityVis]) = tmpResultGUI.(pln.bioParam.quantityVis);
      end
      resultGUI.([pln.bioParam.quantityVis '_' num2str(Cnt,'%d')]) = tmpResultGUI.(pln.bioParam.quantityVis);
      Cnt = Cnt + 1;
   end 
    
end

% remember original fluence weights
resultGUI.w  = w; 




