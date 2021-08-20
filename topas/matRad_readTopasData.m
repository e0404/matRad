function [topasCube,std] = matRad_readTopasData(folder,dij)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if exist('dij')
    calcDoseDirect = false;
else
    calcDoseDirect = true;
end

load([folder filesep 'MCparam.mat']);

%Normalize with histories and particles/weight
correctionFactor = 1e6 / double(MCparam.nbHistoriesTotal); %double(MCparam.nbParticlesTotal) / double(MCparam.nbHistoriesTotal);

cubeDim = MCparam.imageCubeDim;

if calcDoseDirect
    for t = 1:length(MCparam.tallies)
        tname = MCparam.tallies{t};
        topasTally = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        SumVarOverFields = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        
        for f = 1:MCparam.nbFields
            topasSum = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
            topasMeanDiff = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
            for k = 1:MCparam.nbRuns
                genFileName = sprintf('score_%s_field%d_run%d_%s',MCparam.simLabel,f,k,tname);
                switch MCparam.outputType
                    case 'csv'
                        genFullFile = fullfile(folder,[genFileName '.csv']);
                        data{k} = readCsvData(genFullFile,cubeDim);
                    case 'binary'
                        genFullFile = fullfile(folder,[genFileName '.bin']);
                        %                         if iscell(MCparam.scoreReportQuantity)
                        %                             [data{k},topasStd] = readBinData(genFullFile,cubeDim,numel(MCparam.scoreReportQuantity));
                        %                         else
                        data{k} = readBinData(genFullFile,cubeDim);
                        %                         end
                    otherwise
                        error('Not implemented!');
                end
                topasSum = topasSum + data{k};
            end
            
            if strcmp(tname,'physicalDose')
                topasSum = correctionFactor.*topasSum;
                
                % Calculate Standard Deviation from batches
                for k = 1:MCparam.nbRuns
                    topasMeanDiff = topasMeanDiff + (data{k} - topasSum./MCparam.nbRuns).^2;
                end
                % variance of the mean
                topasVarMean = topasMeanDiff./(MCparam.nbRuns - 1)./MCparam.nbRuns;
                % std of the MEAN!
                topasStdMean = sqrt(topasVarMean);
                % std of the SUM
                topasStdSum = topasStdMean * MCparam.nbRuns;
                topasVarSum = topasStdSum.^2;
                
                topasCube.([tname '_std_beam' num2str(f)]) = topasStdSum;
                
                SumVarOverFields = SumVarOverFields + topasVarSum;
            elseif strcmp(tname,'alpha') || strcmp(tname,'beta') || strcmp(tname,'RBE') || strcmp(tname,'LET')
                topasSum = topasSum./ MCparam.nbRuns;
            end
            
            % Accumulate over the fields
            topasTally = topasTally + topasSum;
            % Tally per field
            topasCube.([tname '_beam' num2str(f)]) = topasSum;
            
        end
        if strcmp(tname,'physicalDose')
            topasCube.physicalDose_std = sqrt(SumVarOverFields);
        end
        if ~(strcmp(tname,'alpha') || strcmp(tname,'beta') || strcmp(tname,'RBE'))
            topasCube.(tname) = topasTally;
        end
    end
    
else % if topas dij calculation
    for t = 1:length(MCparam.tallies)
        tname = MCparam.tallies{t};
        topasTally = zeros(cubeDim(1),cubeDim(2),cubeDim(3));
        
        for f = 1:MCparam.nbFields
            topasSum = cell(1,dij.totalNumOfBixels);
            topasSum(1,:) = {zeros(cubeDim(1),cubeDim(2),cubeDim(3))};
            for i = 1:dij.totalNumOfBixels
                for k = 1:MCparam.nbRuns
                    
                    genFileName = sprintf('score_%s_field%d_run%d_%s-ray%i_bixel%i',MCparam.simLabel,f,k,tname,dij.rayNum(i),dij.bixelNum(i));
                    switch MCparam.outputType
                        case 'csv'
                            genFullFile = fullfile(folder,[genFileName '.csv']);
                            data = readCsvData(genFullFile,cubeDim);
                        case 'binary'
                            genFullFile = fullfile(folder,[genFileName '.bin']);
                            data = readBinData(genFullFile,cubeDim);
                        otherwise
                            error('Not implemented!');
                    end
                    topasSum{i} = topasSum{i} + data;
                end
                
                topasSum{i} = correctionFactor.*topasSum{i};
                
                % Tally per field
                topasCube.([tname '_beam' num2str(f)]){i} = topasSum{i};
            end
        end
    end
end



end

function data = readCsvData(csvFile,cubeDim)
data = zeros(cubeDim(2),cubeDim(1),cubeDim(3));
fID = fopen(csvFile,'r');
dataCsv = textscan(fID,'%d %d %d %f','Delimiter',',','CommentStyle','#','CollectOutput',true);
fclose(fID);
ix = sub2ind([cubeDim(1) cubeDim(2) cubeDim(3)],dataCsv{1}(:,2)+1,dataCsv{1}(:,1)+1,dataCsv{1}(:,3)+1);
data(ix) = dataCsv{2};
end

function data = readBinData(binFile,cubeDim)
fID = fopen(binFile);
numOfQuantities = 1;
dataRead = fread(fID,[numOfQuantities,prod(cubeDim)],'double');
data = dataRead(1,:);
% if numOfQuantities ~= 1
%     std = dataRead(2,:);
%     std = reshape(std,cubeDim(2),cubeDim(1),cubeDim(3));
%     std = permute(std,[2 1 3]);
% end
fclose(fID);
data = reshape(data,cubeDim(2),cubeDim(1),cubeDim(3));
data = permute(data,[2 1 3]);
end

