function resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,w)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dose calculation wrapper bypassing dij calculation
% 
% call
%   dij = matRad_calcDoseDirect(ct,stf,pln,cst)
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
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,calcDoseDirect);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst,5000,4,calcDoseDirect);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect);
end

% remember bixel weight
counter = 0;
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
if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
    && strcmp(pln.radiationMode,'carbon')

    ix = resultGUI.physicalDose>0;

    resultGUI.effect     = zeros(ct.cubeDim);
    resultGUI.effect(ix) = dij.mAlphaDose{1}(ix,1) + dij.mSqrtBetaDose{1}(ix,1).^2;

    a_x = zeros(size(resultGUI.physicalDose));
    b_x = zeros(size(resultGUI.physicalDose));

    for i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,4}{1}) = cst{i,5}.alphaX;
            b_x(cst{i,4}{1}) = cst{i,5}.betaX;
        end
    end

    resultGUI.RBExDose = zeros(ct.cubeDim);
    resultGUI.RBExDose(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.effect(ix)) - a_x(ix))./(2.*b_x(ix)));

    resultGUI.alpha    = zeros(ct.cubeDim);
    resultGUI.alpha(ix) = dij.mAlphaDose{1}(ix,1)./resultGUI.physicalDose(ix);

    resultGUI.beta     = zeros(ct.cubeDim);
    resultGUI.beta(ix) = (dij.mSqrtBetaDose{1}(ix,1)./resultGUI.physicalDose(ix)).^2;
    
end

