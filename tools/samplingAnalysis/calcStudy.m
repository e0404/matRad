function calcStudy(examineStructures, multScen, param)
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
[mRealizations,stats, cst, pln, resultCubes,nominalScenario]  = matRad_sampling(ct,stf,cst,pln,resultGUI.w,examineStructures, multScen, param);

%% perform analysis
[structureStat, doseStat] = samplingAnalysis(ct,cst,pln.multScen.subIx,mRealizations,pln.multScen.scenProb);

%% save
filename = 'resultSampling';
save(filename);

%% generate report

cd(param.outputPath)
mkdir(fullfile('report','data'));
mkdir(fullfile('report','data','figures'));
param.outputPath = fullfile('report','data');
copyfile(fullfile(matRadPath,'tools','samplingAnalysis','main_template.tex'),'report/main.tex');

% generate actual latex report
latexReport(ct, cst, pln, nominalScenario, structureStat, param);

cd('report');
executeLatex = 'xelatex main.tex';
system(executeLatex);
system(executeLatex);
system('main.pdf');
