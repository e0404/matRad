function [cstStat, doseStat, param] = matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty sampling analysis function
% 
% call
%   [structureStat, doseStat] = samplingAnalysis(ct,cst,subIx,mSampDose,w)
%
% input
%   ct:             ct cube
%   cst:            matRad cst struct
%   pln:            matRad's pln struct
%   caSampRes       cell array of sampling results depicting plan parameter
%   mSampDose       matrix holding the sampled doses, each row corresponds to
%                   one dose sample
%
% output
%   cstStat         structure-wise statistics (mean, max, percentiles, ...)
%   doseStat        dose-wise statistics (mean, max, percentiles, ...)
%   param           contains additional information about sampling analysis
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

if exist('param','var')
   if ~isfield(param,'logLevel')
       param.logLevel = 1;
   end
   if ~isfield(param,'calcDoseDirect')
      param.calcDoseDirect = false;
   end
else
   param.calcDoseDirect = false;
   param.subIx          = [];
   param.logLevel       = 1;
end

%% check integrity of statistics
param.sufficientStatistics = matRad_checkSampIntegrity(pln.multScen, param);

ix    = pln.subIx;
vProb = pln.multScen.scenProb;

%% generate structurewise statistics
cstStat = struct();
% reassing dvh to stats structure
for i = 1:size(resultGUInomScen.cst,1)
    cstStat(i).name = caSampRes(1).qi(i).name;
    for l = 1:numel(caSampRes)
        % check if structures still match
        if any(~strcmp(cstStat(i).name, {caSampRes(l).dvh(i).name, caSampRes(l).qi(i).name}))
            matRad_dispToConsole('matRad: Error, wrong structure.' , param, 'error');
        end
        cstStat(i).dvh(l).doseGrid     = caSampRes(l).dvh(i).doseGrid;
        cstStat(i).dvh(l).volumePoints = caSampRes(l).dvh(i).volumePoints;
        cstStat(i).qi(l)               = caSampRes(l).qi(i);
        cstStat(i).w(l)                = vProb(l)';
    end  
end
%% calculate mean and std cube
% compute doseMatrix with columns correspond to scenarios


doseStat.meanCube              = zeros(ct.cubeDim);
doseStat.stdCube               = zeros(ct.cubeDim);

doseStat.meanCubeW             = zeros(ct.cubeDim);
doseStat.stdCubeW              = zeros(ct.cubeDim);

doseStat.meanCube(ix)  = mean(mSampDose,2);   
doseStat.stdCube(ix)   = std(mSampDose,1,2);  
doseStat.meanCubeW(ix) = (sum(mSampDose * diag(vProb),2));
doseStat.stdCubeW(ix)  = std(mSampDose,vProb,2);

% gamma cube
doseCube = resultGUInomScen.(pln.bioParam.quantityVis);
if strncmp(pln.bioParam.quantityVis,'RBExD', 5)
    doseStat.gammaAnalysis.cube1Name = 'resultGUInomScen.RBExD';
else
    doseStat.gammaAnalysis.cube1Name = 'resultGUInomScen.physicalDose';
end

doseStat.gammaAnalysis.cube1 = doseCube;
doseStat.gammaAnalysis.cube2 = doseStat.meanCubeW;
doseStat.gammaAnalysis.cube2Name = 'doseStat.meanCubeW';
if ~isfield(param, 'criteria')
    param.criteria = [2 2];
end
matRad_dispToConsole(['matRad: Performing gamma index analysis with parameters', num2str(param.criteria), '[% mm] \n'],param,'info');
doseStat.gammaAnalysis.doseAgreement = param.criteria(1);
doseStat.gammaAnalysis.distAgreement = param.criteria(2);

doseStat.gammaAnalysis.gammaCube = matRad_gammaIndex(doseCube,doseStat.meanCubeW,[ct.resolution.x ct.resolution.y ct.resolution.z],param.criteria);

%% percentiles
if ~isfield(param, 'percentiles')
    param.percentiles = [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99];
end
percentiles     = param.percentiles;
percentileNames = cell(numel(percentiles),1);
% create fieldnames
for i = 1:numel(percentiles)
    percentileNames{i} = ['P',num2str(percentiles(i)*100)];
end
% create table rownames
metric = vertcat({'mean';'min';'max';'std'},percentileNames{:});

% create statstics where structure based results (QI and DVH) are available
for i = 1:size(cst,1)
    cstStat(i).percentiles = percentiles;
    cstStat(i).metric      = metric;
    cstStat(i).dvhStat     = calcDVHStat(cstStat(i).dvh,cstStat(i).percentiles,cstStat(i).w);
    cstStat(i).qiStat      = calcQiStat(cstStat(i).qi,cstStat(i).percentiles,cstStat(i).w);
end


% dvh statistics
    function [dvhStat, doseGrid] = calcDVHStat(dvh,percentiles,w)
        doseGrid = dvh(1).doseGrid;
        dvhMat = NaN * ones(numel(dvh),numel(dvh(1).volumePoints));
        for j = 1:numel(dvh)
            dvhMat(j,:) = dvh(j).volumePoints;            
        end
        % for statistical reasons, treat NaN as 0
        dvhMat(isnan(dvhMat)) = 0;
        
        dvhStat.mean.doseGrid     = doseGrid;
        dvhStat.mean.volumePoints = wMean(dvhMat,w);
        % [~,argmin] = min(dvhStat.mean.volumePoints);
        % dvhStat.mean(2,argmin + 1:end) = NaN;
        
        dvhStat.min.doseGrid     = doseGrid;
        dvhStat.min.volumePoints = min(dvhMat);
        % [~,argmin] = min(dvhStat.min.volumePoints);
        % dvhStat.min(2,argmin + 1:end) = NaN;
        
        dvhStat.max.doseGrid     = doseGrid;
        dvhStat.max.volumePoints = max(dvhMat);
        % [~,argmin] = min(dvhStat.max.volumePoints);
        % dvhStat.max(2,argmin + 1:end) = NaN;
        
        dvhStat.std.doseGrid     = doseGrid;
        dvhStat.std.volumePoints = std(dvhMat,w);

        dvhStat.percDVH = NaN * ones(numel(percentiles),numel(doseGrid));
        
        for j = 1:size(dvhMat,2)
            wQ =  matRad_weightedQuantile(dvhMat(:,j), percentiles, w', false, 'none');
            dvhStat.percDVH(:,j) = wQ;
        end

    end % eof calcDVHStat

    % qi statistics
    function qiStat = calcQiStat(qi,percentiles,w)
        fields = fieldnames(qi);
        % remove name field
        if sum(strcmp('VOIname', fields)) >= 1
            qi = rmfield(qi, 'VOIname');
        end
        fields = fieldnames(qi);
        qiStruct = qi;
        
        % create helper matlab structure which will be converted to table
        qiStatH = struct();
        for j = 1:numel(fields)
            if numel([qiStruct(:).(fields{j})]) == numel(w)
                qiStatH(1).(fields{j}) = wMean([qiStruct(:).(fields{j})],w);
                qiStatH(2).(fields{j}) = min([qiStruct(:).(fields{j})]);
                qiStatH(3).(fields{j}) = max([qiStruct(:).(fields{j})]);
                qiStatH(4).(fields{j}) = std([qiStruct(:).(fields{j})],w);
                wQ = matRad_weightedQuantile([qiStruct(:).(fields{j})], percentiles, w', false, 'none');
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

    function S = wMean(X,w)
        if exist('w','var') || ~isempty(w)
            if isvector(X) && isvector(w)
                S = reshape(w,1,[]) * reshape(X,[],1) / sum(w);
            else
                % row-wise
                S = reshape(w,1,[]) * X ./ sum(w);        
            end

        else
            S = mean(X);
        end
    end

    % check integrity of scenario analysis (i.e. check number of scenarios)
    function statCheck = matRad_checkSampIntegrity(multScen, param)
        if exist('param','var') && isfield(param, 'sufficientStatistics') && ~param.sufficientStatistics
            statCheck = false;
            return
        end
        
        if multScen.totNumScen > 20
            totalNum = true;
        else
            totalNum = false;
        end
        
        statCheck = totalNum; % * .... *
        
    end

end
