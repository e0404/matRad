function resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,w,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dose calculation wrapper bypassing dij calculation
% 
% call
%   resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   w:          optional (if no weights available in stf): bixel weight
%               vector
%   param:      (optional) structure defining additional parameter
%               e.g. param.logLevel
% 
% output
%   resultGUI:  matRad result struct
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define default log level
param.calcDoseDirect = true;
% define certain indices set that should be used for dose calculation
%param.subIx = [];


% copy bixel weight vector into stf struct
if exist('w','var')
    if sum([stf.totalNumOfBixels]) ~= numel(w)
        matRad_dispToConsole('weighting does not match steering information','error')
    end
    counter = 0;
    for i = 1:pln.numOfBeams
        for j = 1:stf(i).numOfRays
            for k = 1:stf(i).numOfBixelsPerRay(j)
                counter = counter + 1;
                stf(i).ray(j).weight(k) = w(counter);
            end
        end
    end
end

if ~isfield(pln,'bioParam')
   % retrieve model parameters
   pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);
end

if ~isfield(pln,'multScen')
   % set plan uncertainties
   [pln] = matRad_setPlanUncertainties(ct,pln);
end

% dose calculation
if strcmp(pln.radiationMode,'photons')
  %dij = matRad_calcPhotonDose(ct,stf,pln,cst,calcDoseDirect);
  dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,5000,4,param);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
  dij = matRad_calcParticleDose(ct,stf,pln,cst,param);
end

% calculate cubes; use uniform weights here, weighting with actual fluence 
% already performed in dij construction 
resultGUI    = matRad_calcCubes(ones(pln.numOfBeams,1),dij,cst);

% remember original fluence weights
resultGUI.w  = w; 

