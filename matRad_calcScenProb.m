function scenProb = matRad_calcScenProb(mu,sigma,samplePos,calcType,probDist,vRange)

if isequal(probDist,'normDist')
    scenProb = 1;
    
    if isequal(calcType,'probBins')
        
        for i = 1:length(mu)
            samplePosSorted = sort(unique(samplePos(i,:)));
            if numel(samplePosSorted) == 1
                continue;
            end
            binWidth        = (samplePosSorted(2) - samplePosSorted(1));
            lowerBinLevel   = samplePos(i,:) - 0.5*binWidth;
            upperBinLevel   = samplePos(i,:) + 0.5*binWidth;
               
            scenProb        = scenProb.*0.5.*(erf((upperBinLevel-mu(i))/(sqrt(2)*sigma(i)))-erf((lowerBinLevel-mu(i))/(sqrt(2)*sigma(i))));
        end
        
    elseif isequal(calcType,'pointwise')
        for i = 1:length(mu)
            scenProb = scenProb .* (1/sqrt(2*pi*sigma(i)^2)*exp(-((samplePos(i,:)-mu(i)).^2./(2*sigma(i)^2))));
        end
    end
    
    % normalize probabilities since we use only a subset of
    % the 3D grid 
    %scenProb = scenProb./sum(scenProb);
    
else
    error('Until now, only normally distributed scenarios implemented')
end
