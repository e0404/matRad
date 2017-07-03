function [resultCubes, percStats] = samplingAnalysis(mRealizations, stats, ct, cst, pln)

% weights
weights = pln.multScen.scenProb;

%% calculate mean and std cube
resultCubes.meanCube              = zeros(ct.cubeDim);
resultCubes.stdCube               = zeros(ct.cubeDim);

resultCubes.meanCubeWeighted      = zeros(ct.cubeDim);
resultCubes.stdCubeWeighted       = zeros(ct.cubeDim);

resultCubes.meanCube(pln.multScen.subIx)         = mean(mRealizations,2);   
resultCubes.stdCube(pln.multScen.subIx)          = std(mRealizations,1,2);  
resultCubes.meanCubeWeighted(pln.multScen.subIx) = (sum(mRealizations * diag(pln.multScen.scenProb),2) )/pln.numOfSamples;
resultCubes.stdCube(pln.multScen.subIx)          = std(mRealizations,pln.multScen.scenProb,2);

%% percentiles
percentiles = [0.1 0.2 0.3 0.7 0.8 0.9];

%percStats = struct(size(cst,1));
for i = 1:size(cst,1)
    percStats(i).percentiles = percentiles;
    percStats(i).percOverPoints = NaN * ones(numel(percentiles),size(stats{1,2}{1,1},2));
end

% dose percentiles
% QI percentiles

% dvh percentiles
dvh = vertcat(stats{:,2});

for j = 1:size(dvh,2)
    singleStructDvh = vertcat(dvh{:,j});
    for k = 1:size(singleStructDvh,2)
        values = singleStructDvh(:,k);
        wQ = weightedQuantile(values, percentiles, weights, false, true);
        percStats(j).percOverPoints(:,k) = wQ';
    end
end

end

