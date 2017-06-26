function resultCubes = samplingAnalysis(mRealizations, ct, pln)

%% calculate mean and std cube
resultCubes.meanCube              = zeros(ct.cubeDim);
resultCubes.stdCube               = zeros(ct.cubeDim);

resultCubes.meanCubeWeighted      = zeros(ct.cubeDim);
resultCubes.stdCubeWeighted       = zeros(ct.cubeDim);

resultCubes.meanCube(pln.multScen.subIx)         = mean(mRealizations,2);   
resultCubes.stdCube(pln.multScen.subIx)          = std(mRealizations,1,2);  
resultCubes.meanCubeWeighted(pln.multScen.subIx) = (sum(mRealizations * diag(pln.multScen.scenProb),2) )/pln.numOfSamples;
resultCubes.stdCube(pln.multScen.subIx)          = std(mRealizations,pln.multScen.scenProb,2);

%% sort cube
[cubesSorted, sortIndex] = sort(mRealizations, 2);

end

