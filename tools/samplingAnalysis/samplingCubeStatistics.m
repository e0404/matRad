function [resultCubes] = samplingCubeStatistics(pln, plnTot, mRealizations, )

resultCubes.meanCube              = zeros(ct.cubeDim);
resultCubes.stdCube               = zeros(ct.cubeDim);

resultCubes.meanCubeWeighted      = zeros(ct.cubeDim);
resultCubes.stdCubeWeighted       = zeros(ct.cubeDim);

resultCubes.meanCube(param.subIx)         = mean(mRealizations,2);   
resultCubes.stdCube(param.subIx)          = std(mRealizations,1,2);  
resultCubes.meanCubeWeighted(param.subIx) = (sum(mRealizations * diag(plnTot.multScen.scenProb),2) )/pln.numOfSamples;
resultCubes.stdCube(param.subIx)          = std(mRealizations,plnTot.multScen.scenProb,2); 


end

