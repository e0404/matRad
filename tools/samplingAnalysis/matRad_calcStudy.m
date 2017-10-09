function matRad_calcStudy(examineStructures, multScen, param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty study wrapper
% 
% call
%   calcStudy(examineStructures, multScen, param)
%
% input
%   examineStructures:  structures which should be examined (can be empty, 
%                       to examine all structures) cube
%   multScen:           parameterset of uncertainty analysis
%   param:              structure defining additional parameter
%                       outputPath
% output
%   (binary)            all results are saved; a pdf report will be generated 
%                       and saved
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
       param.logLevel = 4;
    end   
else
   param.logLevel     = 4;
end

% require minimum number of scenarios to ensure proper statistics
if multScen.numOfRangeShiftScen + sum(multScen.numOfShiftScen) < 12
    warning('You use a very low number of scenarios. Proceeding is not recommended.');
    pause(10);
end

%% load DICOM imported patient
listOfMat = dir('*.mat');
if numel(listOfMat) == 1
  load(listOfMat.name);
else
   matRad_dispToConsole('Ambigous set of .mat files in the current folder (i.e. more than one possible patient or already results available).\n',param,'error');
   return
end

% matRad path
matRadPath = which('matRad.m');
if isempty(matRadPath) 
    matRad_dispToConsole('Please include matRad in your searchpath.',param,'error');
else
    matRadPath = matRadPath(1:(end-8));
end
addpath(fullfile(matRadPath,'tools','samplingAnalysis'));

pln.robOpt = false;
pln.sampling = true;

%% perform calculation
[mRealizations, cst, pln, nominalScenario]  = matRad_sampling(ct,stf,cst,pln,resultGUI.w,examineStructures, multScen, param);

%% perform analysis
[structureStat, doseStat] = matRad_samplingAnalysis(ct,cst,pln.multScen.subIx,mRealizations,pln.multScen.scenProb);

%% save
param.reportPath = fullfile('report','data');
filename = 'resultSampling';
save(filename);

%% generate report

cd(param.outputPath)
mkdir(fullfile('report','data'));
mkdir(fullfile('report','data','figures'));
copyfile(fullfile(matRadPath,'tools','samplingAnalysis','main_template.tex'),fullfile('report','main.tex'));

% generate actual latex report
matRad_latexReport(ct, cst, pln, nominalScenario, structureStat, doseStat, resultGUI, param);

cd('report');
if ispc
    executeLatex = 'xelatex --shell-escape --interaction=nonstopmode main.tex';
elseif isunix
    executeLatex = '/Library/TeX/texbin/xelatex --shell-escape --interaction=nonstopmode main.tex';
end

response = system(executeLatex);
if response == 127 % means not found
    warning('Could not find tex distribution. Please compile manually.');
else
    system(executeLatex);
    system('main.pdf');
end
