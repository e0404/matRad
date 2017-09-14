function resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,w)
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

calcDoseDirect = true;

% copy bixel weight vector into stf struct
if exist('w','var')
    if sum([stf.totalNumOfBixels]) ~= numel(w)
        error('weighting does not match steering information')
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

% dose calculation
dij = cell(length(stf),1);
totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
for i = 1:length(stf)
    if strcmp(pln.radiationMode,'photons')
        dij{i} = matRad_calcPhotonDose(ct,stf(i),pln,cst,calcDoseDirect);
        %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,5000,4,calcDoseDirect);
    elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
        dij{i} = matRad_calcParticleDose(ct,stf(i),pln,cst,calcDoseDirect);
    end
end

% remember bixel weight
counter = 0;
resultGUI.w = NaN * ones(dij{1}.totalNumOfBixels,1);
for i = 1:pln.numOfBeams
    for j = 1:stf(i).numOfRays
        for k = 1:stf(i).numOfBixelsPerRay(j)
            counter = counter + 1;
            resultGUI.w(counter) = stf(i).ray(j).weight(k);
        end
    end
end

for i = 1:(length(stf))
    beamInfo(i).suffix = ['Beam', num2str(i)];
end

% compute phyical dose
sumPhysical = zeros(ct.cubeDim);
for i = 1:length(beamInfo)
    resultGUI.(['physicalDose', beamInfo(i).suffix]) = reshape(full(dij{i}.physicalDose{1}(:,1)),ct.cubeDim);
    sumPhysical = sumPhysical + resultGUI.(['physicalDose', beamInfo(i).suffix]);
end
resultGUI.physicalDose = sumPhysical;

% compute LET if applicable
if isfield(dij{1},'mLETDose')
    sumLET = zeros(ct.cubeDim);
    for i = 1:length(beamInfo)        
        ix = resultGUI.(['physicalDose', beamInfo(i).suffix]) > 0;
        resultGUI.(['LET', beamInfo(i).suffix]) = zeros(dij.dimensions);
        resultGUI.(['LET', beamInfo(i).suffix])(ix) = dij{i}.mLETDose{1}(ix,1)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
        sumLET = sumLET + dij{i}.mLETDose{1}(ix,1);
    end
    ix = resultGUI.physicalDose > 0;
    resultGUI.LET     = zeros(ct.cubeDim);
%   mLETDose = sum(dij{:}.mLETDose{1}(ix,1),4);
    resultGUI.LET(ix) = sumLET./resultGUI.physicalDose(ix);
end
                      
% compute biological cubes
if strcmp(pln.bioOptimization,'const_RBExD')
  for i = 1:length(beamInfo)
       resultGUI.(['RBExDose', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) * dij.RBE;
  end
  resultGUI.RBExDose = resultGUI.physicalDose * dij.RBE;  
  
elseif strcmp(pln.bioOptimization,'LEMIV_effect') || strcmp(pln.bioOptimization,'LEMIV_RBExD')
    
    sumAlpha = dij{1}.mAlphaDose{1}(:,1) * 0;
    sumSqrtBeta = dij{1}.mSqrtBetaDose{1}(:,1) * 0;
    for i = 1:length(beamInfo)
        ix = resultGUI.(['physicalDose', beamInfo(i).suffix]) > 0;

        resultGUI.(['effect', beamInfo(i).suffix]) = zeros(ct.cubeDim);
        resultGUI.(['RBExDose', beamInfo(i).suffix]) = zeros(ct.cubeDim);
        resultGUI.(['alpha', beamInfo(i).suffix]) = zeros(ct.cubeDim);
        resultGUI.(['beta', beamInfo(i).suffix]) = zeros(ct.cubeDim);
        
        resultGUI.(['effect', beamInfo(i).suffix])(ix) = dij{i}.mAlphaDose{1}(ix,1) + dij{i}.mSqrtBetaDose{1}(ix,1).^2;

        a_x = zeros(size(resultGUI.(['physicalDose', beamInfo(i).suffix])));
        b_x = zeros(size(resultGUI.(['physicalDose', beamInfo(i).suffix])));

        for j = 1:size(cst,1)
            % Only take OAR or target VOI.
            if isequal(cst{j,3},'OAR') || isequal(cst{j,3},'TARGET') 
                a_x(cst{j,4}{1}) = cst{j,5}.alphaX;
                b_x(cst{j,4}{1}) = cst{j,5}.betaX;
            end
        end
        
      resultGUI.(['RBExDose', beamInfo(i).suffix])(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) ...
            .* resultGUI.(['effect', beamInfo(i).suffix])(ix)) - a_x(ix))./(2.*b_x(ix)));
      
      resultGUI.(['alpha', beamInfo(i).suffix])(ix) = dij{i}.mAlphaDose{1}(ix,1) ./ resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
      resultGUI.(['beta', beamInfo(i).suffix])(ix) = (dij{i}.mSqrtBetaDose{1}(ix,1) ./ resultGUI.(['physicalDose', beamInfo(i).suffix])(ix)).^2;
      sumAlpha = sumAlpha + dij{i}.mAlphaDose{1}(:,1);
      sumSqrtBeta = sumSqrtBeta + dij{i}.mSqrtBetaDose{1}(:,1);
    end
    
    ix = resultGUI.physicalDose > 0;
    
    resultGUI.effect = zeros(ct.cubeDim);
    resultGUI.RBExDose = zeros(ct.cubeDim);
    resultGUI.alpha = zeros(ct.cubeDim);
    resultGUI.beta = zeros(ct.cubeDim);
    
    resultGUI.effect(ix) = sumAlpha(ix) + sumSqrtBeta(ix).^2;
    
    resultGUI.RBExDose(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.effect(ix)) - a_x(ix))./(2.*b_x(ix)));
    resultGUI.alpha(ix) = sumAlpha(ix)./resultGUI.physicalDose(ix);
    resultGUI.beta(ix) = (sumSqrtBeta(ix)./resultGUI.physicalDose(ix)).^2;
    
end
