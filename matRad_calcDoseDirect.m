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

if exist('param','var')
   if ~isfield(param,'logLevel')
      param.logLevel = 1;
   end
else
   param.logLevel = 1;
end

% copy bixel weight vector into stf struct
if exist('w','var')
    if sum([stf.totalNumOfBixels]) ~= numel(w)
        matRad_dispToConsole('weighting does not match steering information',param,'error')
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

% disable robOpt
pln.robOpt   = false;

% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);

% set plan uncertainties for robust optimization
[pln] = matRad_setPlanUncertainties(ct,pln);

% define certain indices set that should be used for dose calculation
param.subIx = [];

% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,param);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,5000,4,calcDoseDirect);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,param);
end

% remember bixel weight
counter = 0;
resultGUI.w = NaN * ones(dij.totalNumOfBixels,1);
for i = 1:pln.numOfBeams
    for j = 1:stf(i).numOfRays
        for k = 1:stf(i).numOfBixelsPerRay(j)
            counter = counter + 1;
            resultGUI.w(counter) = stf(i).ray(j).weight(k);
        end
    end
end

% compute phyical dose
resultGUI.physicalDose = reshape(full(dij.physicalDose{1}(:,1)),ct.cubeDim);

% compute LET if applicable
if isfield(dij,'mLETDose')
    
    ix = resultGUI.physicalDose>0;
    
    resultGUI.LET     = zeros(ct.cubeDim);
    resultGUI.LET(ix) = dij.mLETDose{1}(ix,1)./resultGUI.physicalDose(ix);
    
end
                      
% compute biological cubes
if isfield(dij,'RBE')

    resultGUI.RBExDose = resultGUI.physicalDose * dij.RBE;
    
elseif isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')

    ix = resultGUI.physicalDose>0;

    resultGUI.effect     = zeros(ct.cubeDim);
    resultGUI.effect(ix) = dij.mAlphaDose{1}(ix,1) + dij.mSqrtBetaDose{1}(ix,1).^2;

    resultGUI.RBExDose     = zeros(ct.cubeDim);
    resultGUI.RBExDose(ix) = ((sqrt(dij.alphaX(ix).^2 + 4 .* dij.betaX(ix) .* resultGUI.effect(ix)) - dij.alphaX(ix))./(2.*dij.betaX(ix)));

    resultGUI.alpha     = zeros(ct.cubeDim);
    resultGUI.alpha(ix) = dij.mAlphaDose{1}(ix,1)./resultGUI.physicalDose(ix);

    resultGUI.beta     = zeros(ct.cubeDim);
    resultGUI.beta(ix) = (dij.mSqrtBetaDose{1}(ix,1)./resultGUI.physicalDose(ix)).^2;
    
end

