function [ Stats ] = matRad_unitTestRun(NameOfUnitTestResult,SaveResultToDisk)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create a reference unit test result
% 
% call
%   Stats = matRad_unitTestRun('unitTest_result_C_001_BOXPHANTOM',0)
%
% input
%   NameOfUnitTestResult:   string determing a unitTest reference result 
%   SaveResultToDisk:       Boolean if the result should be saved in a
%                           seperated text file
%
% output
%   Stats:                  Struct containg a statistical evaulation of the
%                           unit test result
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rootDir = pwd;
rootDir = rootDir(1:end-9);  % /unitTest has exactly nine chars
addpath(rootDir);

%% check how many plans for this patient data already exist 
load(NameOfUnitTestResult);
load(pln.patientname);
RunInfo.Status = 'undefined';
try
    stf = matRad_generateStf(ct,cst,pln);
    if strcmp(pln.radiationMode,'photons')
        dij = matRad_calcPhotonDose(ct,stf,pln,cst,0);
    elseif strcmp(pln.radiationMode,'protons') || ...
            strcmp(pln.radiationMode,'carbon')
        dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
    end
    
    resultGUI = matRad_fluenceOptimization(dij,cst,pln,0);
    
    %% sequencing
    if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
        %resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,5);
        resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,5);
    end

    %% DAO
    if strcmp(pln.radiationMode,'photons') && pln.runDAO
       resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,1);
       matRad_visApertureInfo(resultGUI.apertureInfo);
    end
    %% extract dij from iso center slice according to the sample spacing
    resultGUI.dijIsoCenter = matRad_unitTestGetIsoCenterSlice(pln,dij,pln.dijSpacing);

    RunInfo.Status = 'run successfull';
catch error
    RunInfo.error = error;
    RunInfo.Status = 'run failed';
end


Stats = matRad_verifyResultUnitTest(ExpResultGUI,resultGUI);
StatsOverall.max  = 0;
StatsOverall.min  = 0;
StatsOverall.mean = 0;
StatsOverall.std  = 0;
StatsOverall = matRad_unitTestGetOverallStats(StatsOverall,Stats);

matRad_unitTestDispResult(StatsOverall,Stats,NameOfUnitTestResult,RunInfo,...
                          SaveResultToDisk);
end




%% displays the result in the command window and saves the detailed statistic structure on disk
function  matRad_unitTestDispResult(StatsOverall,Stats,TestCaseName,RunInfo,saveStats)

%%get overall values

str1 = ('__________________________________________________________________');
str2 = (['testcase ' TestCaseName ': ' RunInfo.Status]);
str3 = ('overall statistics are as follows: ');
str4 = ([' std: ' num2str(StatsOverall.std) ...
         ' min: ' num2str(StatsOverall.min) ...
         ' max: ' num2str(StatsOverall.max) ...
         ' mean: ' num2str(StatsOverall.mean)]);             
str5 = ('__________________________________________________________________');

disp(str1);
disp(str2);
disp(str3);
disp(str4);
disp(str5);


if saveStats
    Name = [TestCaseName '_' datestr(now,'dd-mm-yyyy HH-MM-SS_')  '.mat'];
    save(Name,'Stats');
end

end


%% function extracts dij from voxels belonging to the iso center slice 
function dijIsoCenter = matRad_unitTestGetIsoCenterSlice(pln,dij,Spacing)

StartIdx = pln.voxelDimensions(1)*pln.voxelDimensions(2)*(pln.voxelDimensions(3)-1);
EndIdx   = pln.voxelDimensions(1)*pln.voxelDimensions(2)*(pln.voxelDimensions(3)); 
LinearIdx = StartIdx:Spacing:EndIdx;
dijIsoCenter = dij.physicalDose(LinearIdx,:);

end

%% functions extracts overall statistics from Stats struct to enable a quick review
function StatsOverall = matRad_unitTestGetOverallStats(StatsOverall,Stats)


fNames = fieldnames(Stats);

    for i = 1:numel(fNames)
        for j=1:length(Stats)
            if sum(strcmp((fNames{i,1}),{'min','max','std','mean'}))>0

                if StatsOverall.min < Stats.min
                    StatsOverall.min = Stats.min;
                end
                if StatsOverall.max < Stats.max
                    StatsOverall.max = Stats.max;
                end
                if StatsOverall.std < Stats.std
                    StatsOverall.std = Stats.std;
                end
                if StatsOverall.mean < Stats.mean
                    StatsOverall.mean = Stats.mean;
                end
                break;
            else

                 StatsOverall = matRad_unitTestGetOverallStats(StatsOverall,Stats(j).(fNames{i,1}));

            end
        end
    end

end