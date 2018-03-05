function matRad_calcStudy(structSel,multScen,matPatientPath,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty study wrapper
% 
% call
%   calcStudy(structSel,multScen,param)
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
function [treatmentSimulation, scenContainer, nominalScenario] = mergeSplittedScen(scenContainerParts, nomScenParts)
    % Merge nominal scenario by cumulating on first part. Remember that all 
    % scenario class is subclass of abstract handle class, therefore 
    % nomScenParts{1} changes.
    
    nominalScenario = nomScenParts{1};
    for k = 2:numel(nomScenParts)
        nominalScenario.cumulateDose(nomScenParts{k}.dose);
    end
    scenContainer = scenContainerParts{1};
    for j = 1:numel(scenContainer)
        for k = 2:numel(scenContainerParts)
            scenContainer{j}.cumulateDose(scenContainerParts{k}{j}.dose);
        end
    end
    treatmentSimulation = MatRadSimulation(pln.bioParam.quantityVis, nominalScenario, ct, cst, pln, pln.numOfFractions);
    for j = 1:numel(scenContainer)
        treatmentSimulation.initNewScen(scenContainer{j})
    end
end

function stf = copyWeightsToStf(stf, w)
    if sum([stf.totalNumOfBixels]) ~= numel(w)
        error('weighting does not match steering information')
    end
    counter = 0;
    for ii = 1:numel(stf)
        for j = 1:stf(ii).numOfRays
            for k = 1:stf(ii).numOfBixelsPerRay(j)
                counter = counter + 1;
                stf(ii).ray(j).weight(k) = w(counter);
            end
        end
    end
end

if exist('param','var')
    if ~isfield(param,'logLevel')
       param.logLevel = 4;
    end   
else
   param.logLevel     = 4;
end

if ~isfield(param,'outputPath')
    param.outputPath = mfilename('fullpath');
end
if ~isfield(param, 'percentiles')
    param.percentiles = [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99];
end
if ~isfield(param, 'criteria')
    param.gammaCriteria = [2 2];
end

% require minimum number of scenarios to ensure proper statistics
if multScen.totNumScen < 20 % multScen.numOfRangeShiftScen + sum(multScen.numOfShiftScen) < 20
    matRad_dispToConsole('Detected a low number of scenarios. Proceeding is not recommended.',param,'warning');
    param.sufficientStatistics = false;
    pause(1);
end

%% load DICOM imported patient
cst = [];
ct = [];
pln = [];
dij = [];
resultGUI = [];
stf = [];

if exist('matPatientPath', 'var') && ~isempty(matPatientPath) && exist('matPatientPath','file') == 2
    load(matPatientPath)
else
    listOfMat = dir('*.mat');
    if numel(listOfMat) == 1
      load(listOfMat.name);
    else
       matRad_dispToConsole('Ambigous set of .mat files in the current folder (i.e. more than one possible patient or already results available).\n',param,'error');
       return
    end
end

stf = copyWeightsToStf(stf, resultGUI.w);

% check if nominal workspace is complete
if ~(exist('ct','var') && exist('cst','var') && exist('stf','var') && exist('pln','var') && exist('resultGUI','var'))
    matRad_dispToConsole('Nominal workspace for sampling is incomplete.\n',param,'error');
end

% matRad path
matRadPath = which('matRad.m');
if isempty(matRadPath) 
    matRad_dispToConsole('Please include matRad in your searchpath.',param,'error');
else
    matRadPath = matRadPath(1:(end-8));
end
addpath(fullfile(matRadPath,'tools','samplingAnalysis'));

% calculate RBExDose
if ~isfield(pln, 'bioParam')
    if strcmp(pln.radiationMode, 'protons')
        pln.bioOptimization = 'constRBE_RBExD';
    elseif strcmp(pln.radiationMode, 'carbon')
        pln.bioOptimization = 'LEM_RBExD';
    end
    pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);
end
    

pln.robOpt   = false;
pln.sampling = true;

%% perform calculation and save
multScen = multScen.matRad_createValidInstance();

tic
if strcmp(multScen.deepestCorrelationBreak, 'fraction')
    [treatmentSimulation, scenContainer, pln, resultGUInomScen, nomScen] = matRad_sampling(ct,stf,cst,pln,resultGUI.w,structSel,multScen,param);
elseif strcmp(multScen.deepestCorrelationBreak, 'beam')
    scenContainerParts = cell(numel(stf),1);
    nomScenParts = cell(numel(stf),1);
    for i = 1:numel(stf)
        stfPart = stf(i);
        plnPart = pln;
        plnPart.couchAngles = pln.couchAngles(i);
        plnPart.gantryAngles = pln.gantryAngles(i);
        plnPart.isoCenter = pln.isoCenter(i,:);
        plnPart.numOfBeams = 1;
        if strcmp(multScen.rangeShiftCorrelationBreak, 'beam') && strcmp(multScen.isoShiftCorrelationBreak, 'beam')
            multScen = multScen.matRad_recreateInstance;
        elseif strcmp(multScen.isoShiftCorrelationBreak, 'beam')
            multScen = multScen.matRad_recreateInstance('range');
        elseif strcmp(multScen.rangeShiftCorrelationBreak, 'beam')
            multScen = multScen.matRad_recreateInstance('setup');
        else
            erorr('Incosistency');
        end
        
        [~, scenContainerPart, ~, ~, nomScenPart] = matRad_sampling(ct, stfPart, cst, plnPart, [],structSel, multScen, param);
        nomScenParts{i} = nomScenPart;
        scenContainerParts{i} = scenContainerPart;
    end
    
    [treatmentSimulation, scenContainer, nomScen] = mergeSplittedScen(scenContainerParts, nomScenParts);
end
param.computationTime = toc;

param.reportPath = fullfile('report','data');
filename         = 'resultSampling';
save(filename, '-v7.3');

%% perform analysis 
% start here loading resultSampling.mat if something went wrong during analysis or report generation
treatmentSimulation.runAnalysis(param.gammaCriteria, param.percentiles);

% %% generate report
% listOfQI = {'mean', 'std', 'max', 'min', 'D_2', 'D_5', 'D_50', 'D_95', 'D_98'};
% 
% cd(param.outputPath)
% mkdir(fullfile('report','data'));
% mkdir(fullfile('report','data','frames'));
% mkdir(fullfile('report','data','figures'));
% copyfile(fullfile(matRadPath,'tools','samplingAnalysis','main_template.tex'),fullfile('report','main.tex'));
% 
% % generate actual latex report
% matRad_latexReport(ct, cst, pln, resultGUInomScen, structureStat, doseStat, mSampDose, listOfQI, param);
% 
% cd('report');
% if ispc
%     executeLatex = 'xelatex --shell-escape --interaction=nonstopmode main.tex';
% elseif isunix
%     executeLatex = '/Library/TeX/texbin/xelatex --shell-escape --interaction=nonstopmode main.tex';
% end
% 
% response = system(executeLatex);
% if response == 127 % means not found
%     warning('Could not find tex distribution. Please compile manually.');
% else
%     system(executeLatex);
% end
end
