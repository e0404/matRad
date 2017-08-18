function structureStat = samplingAnalysis(cst,w)

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
percentiles = [0.02 0.1 0.2 0.3 0.7 0.8 0.9 0.98];
percentileNames = cell(numel(percentiles),1);
% create fieldnames
for i = 1:numel(percentiles)
    percentileNames{i} = ['P',num2str(percentiles(i)*100)];
end
% create table rownames
metric = vertcat({'mean';'min';'max';'std'},percentileNames{:});

% create statstics where structure based results (QI and DVH) are available
for i = 1:size(cst,1)
    structureStat(i).name = cst{i,2};
    structureStat(i).dvh = cst{i,8};
    structureStat(i).qi = cst{i,9};
    structureStat(i).percentiles = percentiles;
    structureStat(i).metric = metric;
    
    structureStat(i).dvhStat = calcDVHStat(cst{i,8},percentiles,w);
    structureStat(i).qiStat = calcQiStat(cst{i,9},percentiles,w);
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
        dvhStat.min(1,:) = doseGrid;
        dvhStat.min(2,:) = min(dvhMat);
        dvhStat.max(1,:) = doseGrid;
        dvhStat.max(2,:) = max(dvhMat);
        
        dvhStat.std(1,:) = doseGrid;
        dvhStat.std(2,:) = std(dvhMat,w);

        dvhStat.percDVH = NaN * ones(numel(percentiles),size(dvh{1},2));
        
        for j = 1:size(dvhMat,2)
            wQ =  weightedQuantile(dvhMat(:,j), percentiles, w', false, true);
            dvhStat.percDVH(:,j) = wQ;
        end
    end % eof calcDVHStat

    % qi statistics
    function qiStat = calcQiStat(qi,percentiles,w)
        fields = fieldnames(qi{1});
        % convert to struct (maybe change structure to struct later)
        qiStruct(1) = qi{1};
        for j = 2:numel(qi)
            qiStruct(j) = qi{j};
        end
        
        % create helper matlab structure which will be converted to table
        qiStatH = struct();
        for j = 1:numel(fields)
            if numel([qiStruct(:).(fields{j})]) == numel(w)
                qiStatH(1).(fields{j}) = wMean([qiStruct(:).(fields{j})],w);
                qiStatH(2).(fields{j}) = min([qiStruct(:).(fields{j})]);
                qiStatH(3).(fields{j}) = max([qiStruct(:).(fields{j})]);
                qiStatH(4).(fields{j}) = std([qiStruct(:).(fields{j})],w);
                wQ = weightedQuantile([qiStruct(:).(fields{j})], percentiles, w', false, true);
                for k = 1:numel(wQ)
                    sIx = k + 4;
                    qiStatH(sIx).(fields{j}) = wQ(k);
                end
            else
                for k = 1:(4 + numel(percentiles))
                    qiStatH(k).(fields{j}) = [];
                end
            end
        end
        qiStat = struct2table(qiStatH);
        qiStat.Properties.RowNames = metric;
    end % eof calcQiStat

end

