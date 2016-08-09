function scenProb = matRad_calcScenProb(mu,sigma,samplePosition,calcType,probDist)

if isequal(probDist,'normDist')
    if isequal(calcType,'probBins')
        
        binWidth      = (samplePosition(2) - samplePosition(1));
        lowerBinLevel = samplePosition - 0.5*binWidth;
        upperBinLevel = samplePosition + 0.5*binWidth;
               
        scenProb      = 0.5*(erf((upperBinLevel-mu)/(sqrt(2)*sigma))-erf((lowerBinLevel-mu)/(sqrt(2)*sigma)));
        
    else
        error('Until now, only probBins implemented')
    end
    
else
    error('Until now, only normally distributed scenarios implemented')
end
