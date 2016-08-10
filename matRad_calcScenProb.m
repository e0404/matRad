function scenProb = matRad_calcScenProb(mu,sigma,samplePos,calcType,probDist)

if isequal(probDist,'normDist')
    if isequal(calcType,'probBins')
        
        scenProb = 1;
        
        for i = 1:length(mu)
            samplePosSorted = sort(unique(samplePos(i,:)));
            binWidth        = (samplePosSorted(2) - samplePosSorted(1));
            lowerBinLevel   = samplePos(i,:) - 0.5*binWidth;
            upperBinLevel   = samplePos(i,:) + 0.5*binWidth;
               
            scenProb        = scenProb.*0.5.*(erf((upperBinLevel-mu(i))/(sqrt(2)*sigma(i)))-erf((lowerBinLevel-mu(i))/(sqrt(2)*sigma(i))));
        end
        
        % normalize probabilities since we we use only a small subset of
        % the 3D grid
        scenProb = scenProb./sum(scenProb);
        
    else
        error('Until now, only probBins implemented')
    end
    
else
    error('Until now, only normally distributed scenarios implemented')
end
