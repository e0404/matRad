function scenProb = matRad_calcScenProb(mu,sigma,samplePos,calcType,probDist)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_calcScenProb provides different ways of calculating the probability 
% of occurance of individual scenarios
% 
% call
%   scenProb = matRad_calcScenProb(mu,sigma,samplePos,calcType,probDist)
%
% input
%   mu:             mean of the distrubtion
%   sigma:          standard deviation of the distribution
%   calcType:       can be set to 
%                   (i)  probBins to calculate the accumulated occurance probability in a certain bin width
%                   (ii) pointwise to calculate the pointwise occurance probability
%   probDist:       identifier for the underlying probability distribution
%                   (i) normDist
%
% output
%   scenProb:       occurance probability of the specified scenario
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isequal(probDist,'normDist')
   
    scenProb = 1;
    
    if isequal(calcType,'probBins')
        
        for i = 1:length(mu)
            samplePosSorted = sort(unique(samplePos(:,i)));
            if numel(samplePosSorted) == 1 || sigma(i) == 0
                continue;
            end
            binWidth        = (samplePosSorted(2) - samplePosSorted(1));
            lowerBinLevel   = samplePos(:,i) - 0.5*binWidth;
            upperBinLevel   = samplePos(:,i) + 0.5*binWidth;
               
            scenProb        = scenProb.*0.5.*(erf((upperBinLevel-mu(i))/(sqrt(2)*sigma(i)))-erf((lowerBinLevel-mu(i))/(sqrt(2)*sigma(i))));
        end
        
    elseif isequal(calcType,'pointwise')
        for i = 1:length(mu)
            scenProb = scenProb .* (1/sqrt(2*pi*sigma(i)^2)*exp(-((samplePos(:,i)-mu(i)).^2./(2*sigma(i)^2))));
        end
    end
    
    % normalize probabilities since we use only a subset of
    % the 3D grid 
    scenProb = scenProb./sum(scenProb);

elseif isequal(probDist,'equalProb')
   
   numScen  = size(samplePos,1);
   scenProb = repmat(1/numScen,1,numScen);
   
else
    matRad_dispToConsole('Until now, only normally distributed scenarios implemented','error')
end
