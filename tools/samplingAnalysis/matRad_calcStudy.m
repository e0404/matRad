function matRad_calcStudy(multScen,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty study wrapper
%
% call
%   matRad_calcStudy(structSel,multScen,matPatientPath,param)
%
% input
%   structSel:          structures which should be examined (can be empty,
%                       to examine all structures) cube
%   multScen:           parameterset of uncertainty analysis
%   matPatientPath:     (optional) absolut path to patient mat file. If
%                       empty mat file in current folder will be used
%   param:              structure defining additional parameter
%                       outputPath
% output
%   (binary)            all results are saved; a pdf report will be generated
%                       and saved
%
%
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

matRad_cfg = MatRad_Config.instance();

p = inputParser;
p.addRequired('multScen',@(x) isa(x,'matRad_multScen'));
p.addParameter('SelectStructures',cell(0),@iscellstr);
p.addParameter('OutputPath',mfilename('fullpath'),@isfolder);
p.addParameter('PatientMatFile','',@isfile);
p.addParameter('ListOfQI',{'mean', 'std', 'max', 'min', 'D_2', 'D_5', 'D_50', 'D_95', 'D_98'},@iscellstr);
p.addParameter('OperatorName','matRad User',@(x) isstring(x) || ischar(x));

%

p.parse(multScen,varargin{:});
multScen = p.Results.multScen;
outputPath = p.Results.OutputPath;
structSel = p.Results.SelectStructures;
matPatientPath = p.Results.PatientMatFile;
listOfQI = p.Results.ListOfQI;
operator = p.Results.OperatorName;


% require minimum number of scenarios to ensure proper statistics
if multScen.numOfRangeShiftScen + sum(multScen.numOfShiftScen) < 20
    matRad_cfg.dispWarning('Detected a low number of scenarios. Proceeding is not recommended.');
    sufficientStatistics = false;
    pause(1);
else
    sufficientStatistics = true;
end

%% load DICOM imported patient or run from workspace
if exist('matPatientPath', 'var') && ~isempty(matPatientPath) && exist('matPatientPath','file') == 2
    load(matPatientPath);
else
    try
        ct          = evalin('base','ct');
        cst         = evalin('base','cst');
        stf         = evalin('base','stf');
        pln         = evalin('base','pln');
        resultGUI   = evalin('base','resultGUI');
    catch
        matRad_cfg.dispError('Workspace for sampling is incomplete.');
    end
end

% check if nominal workspace is complete
if ~(exist('ct','var') && exist('cst','var') && exist('stf','var') && exist('pln','var') && exist('resultGUI','var'))
    matRad_cfg.dispError('Workspace for sampling is incomplete.');
end

% calculate RBExDose
if ~isfield(pln, 'bioParam')
    if strcmp(pln.radiationMode, 'protons')
        pln.bioOptimization = 'RBExD';
        pln.model = 'constRBE';
    elseif strcmp(pln.radiationMode, 'carbon')
        pln.bioOptimization = 'RBExD';
        pln.model = 'LEM';
    end
    pln.bioParam = matRad_bioModel(pln.radiationMode, pln.bioOptimization, pln.model);
end


pln.robOpt   = false;
pln.sampling = true;

%% perform calculation and save
tic
[caSampRes, mSampDose, pln, resultGUInomScen]  = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,multScen);
computationTime = toc;

filename         = 'resultSampling';
save(filename, '-v7.3');

%% perform analysis
% start here loading resultSampling.mat if something went wrong during analysis or report generation
[structureStat, doseStat, meta] = matRad_samplingAnalysis(ct,cst,pln,caSampRes,mSampDose,resultGUInomScen);

%% generate report
matRadPath = matRad_cfg.matRadRoot;
reportPath = 'report';
   
mkdir([outputPath filesep reportPath]);
    
copyfile(fullfile(matRadPath,'tools','samplingAnalysis','main_template.tex'),fullfile(outputPath,reportPath,'main.tex'));
    
% generate actual latex report
success = matRad_latexReport([outputPath filesep reportPath],ct, cst, pln, resultGUInomScen, structureStat, doseStat, mSampDose, listOfQI,...
    'ComputationTime',computationTime,...
    'SufficientStatistics',sufficientStatistics,...
    'OperatorName',operator);


if success
    open(fullfile([outputPath filesep reportPath],'main.pdf'));
    
else
     matRad_cfg.dispError('Report PDF can not be opened...');
end


   
    

