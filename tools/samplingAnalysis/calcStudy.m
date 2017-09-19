function calcStudy(examineStructures, multScen, param)
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

%% load DICOM imported patient
listOfMat = dir('*.mat');
if numel(listOfMat) == 1
  load(listOfMat.name);
else
   matRad_dispToConsole('Ambigous set of .mat files in the current folder (i.e. more than one possible patient).',param,'error');
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
[structureStat, doseStat] = samplingAnalysis(ct,cst,pln.multScen.subIx,mRealizations,pln.multScen.scenProb);

%% save
filename = 'resultSampling';
save(filename);

%% generate report

cd(param.outputPath)
mkdir(fullfile('report','data'));
mkdir(fullfile('report','data','figures'));
param.reportPath = fullfile('report','data');
copyfile(fullfile(matRadPath,'tools','samplingAnalysis','main_template.tex'),'report/main.tex');

% generate actual latex report
latexReport(ct, cst, pln, nominalScenario, structureStat, resultGUI, param);

cd('report');
executeLatex = 'xelatex -shell-escape main.tex';
system(executeLatex);
system(executeLatex);
system('main.pdf');
