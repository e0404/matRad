function [structureStat] = samplingAnalysis(cst,w)%(mRealizations, stats, ct, cst, pln)

% weights
% weights = pln.multScen.scenProb;
% 
% % calculate mean and std cube
% resultCubes.meanCube              = zeros(ct.cubeDim);
% resultCubes.stdCube               = zeros(ct.cubeDim);
% 
% resultCubes.meanCubeWeighted      = zeros(ct.cubeDim);
% resultCubes.stdCubeWeighted       = zeros(ct.cubeDim);
% 
% resultCubes.meanCube(pln.multScen.subIx)         = mean(mRealizations,2);   
% resultCubes.stdCube(pln.multScen.subIx)          = std(mRealizations,1,2);  
% resultCubes.meanCubeWeighted(pln.multScen.subIx) = (sum(mRealizations * diag(pln.multScen.scenProb),2) )/pln.numOfSamples;
% resultCubes.stdCube(pln.multScen.subIx)          = std(mRealizations,pln.multScen.scenProb,2);

%% percentiles
percentiles = [0.1 0.2 0.3 0.7 0.8 0.9];
% create fieldnames 

% create statstics where structure based results (QI and DVH) are available
for i = 1:size(cst,1)
    structureStat(i).name = cst{i,2};
    structureStat(i).dvh = cst{i,8};
    structureStat(i).qi = cst{i,9};
    structureStat(i).percentiles = percentiles;
    
    structureStat(i).dvhStat = calcDVHStat(cst{i,8},percentiles,w);
end


% dvh statistics

    function dvhStat = calcDVHStat(dvh,percentiles,w)
        doseGrid = dvh{1}(1,:);
        dvhMat = NaN * ones(numel(dvh),size(dvh{1},2));
        for j = 1:numel(dvh)
            dvhMat(j,:) = dvh{j}(2,:);
        end
        dvhStat.mean(1,:) = doseGrid;
        dvhStat.mean(2,:) = wMean(dvhMat,w);

        dvhStat.percDVH = NaN * ones(numel(percentiles),size(dvh{1},2));
        
        for j = 1:size(dvhMat,2)
            wQ =  weightedQuantile(dvhMat(:,j), percentiles, w', false, true);
            dvhStat.percDVH(:,j) = wQ;
        end
        
        dvhStat.std(1,:) = doseGrid;
        dvhStat.std(2,:) = std(dvhMat,w);
        
    end

% qi statistics

    function qiStat = calcQiStat(qi,w)
        fields = fieldnames(qi{1});
        % convert to struct (maybe change structure to struct later)
        qiStruct(1) = qi{1};
        for j = 2:numel(qi)
            qiStruct(j) = qi{j};
        end
    end

end

