function [treatmentSimulation, scenContainer, pln, resultGUInomScen, nomScen]  = matRad_sampling(ct,stf,cst,pln,w,structSel,multScen,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_randomSampling enables sampling multiple treatment scenarios
% 
% call
%   [cst,pln] = matRad_setPlanUncertainties(ct,cst,pln)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   w:          optional (if no weights available in stf): bixel weight
%               vector
%   param:      (optional) structure defining additional parameter
%               e.g. param.logLevel
%
% output
%   caSampRes:         cell array of sampling results depicting plan parameter
%   mSampDose:         matrix holding the sampled doses, each row corresponds to
%                      one dose sample
%   pln:               matRad pln struct containing sampling information
%   resultGUInomScen:  resultGUI struct of the nominal scenario
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
    % default: dose influence matrix computation
   if ~isfield(param,'calcDoseDirect')
      param.calcDoseDirect = false;
   end
    
else
   param.calcDoseDirect = false;
   param.subIx          = [];
   param.logLevel       = 1;
end


% save nonSampling pln for nominal scenario calculation and add dummy fields
plnNominal = pln;
% create nominal scenario
plnNominal.multScen = matRad_multScen(ct,'nomScen'); 

% either use existing multScen struct or create new one
if exist('multScen','var') && ~isempty(multScen)
    pln.multScen = multScen;
else
    % create random scenarios for sampling
    pln.multScen = matRad_multScen(ct,'rndScen'); % 'impSamp' or 'wcSamp'
end

matRad_dispToConsole(['Using ' num2str(pln.multScen.totNumScen) 'samples in total \n'],param,'info')

V = [];
% define voxels for sampling
if ~exist('structSel', 'var') || sum(size(structSel)) == 0
    V = [cst{:,4}];
else
    for i=1:size(cst,1)
        for j = 1:numel(structSel)
            if strcmp(structSel{j}, cst{i,2})
                V = [V {cst{i,4}{1}}];
            end
        end
    end
end

% final voxel subset for sampling
param.subIx = unique(vertcat(V{:}));

% disable structures for DVH plotting which are not completely in subIx
for i = 1:size(cst,1)
    if ~all(ismember(cst{i,4}{1}, param.subIx))
        cst{i,5}.Visible = false;
    end
end

% define variable for storing scenario doses
mSampDose   = single(zeros(numel(param.subIx),pln.multScen.totNumScen,1));
StorageInfo = whos('mSampDose');
matRad_dispToConsole(['matRad: Realizations variable will need: ' num2str(StorageInfo.bytes/1e9) ' GB \n'],param,'info');

% check if parallel toolbox is installed and license can be checked out
try
    [FlagParallToolBoxLicensed,msg]  = license('checkout','Distrib_Computing_Toolbox');
    if ~FlagParallToolBoxLicensed
        matRad_dispToConsole(['Could not check out parallel computing toolbox. \n'],param,'warning');
    end
    
catch
    FlagParallToolBoxLicensed  = false;
end

%% calculate nominal scenario
nomScenTimer     = tic;
resultGUInomScen = matRad_calcDoseDirect(ct,stf,plnNominal,cst,w,param);
nomScenTime      = toc(nomScenTimer);
matRad_dispToConsole(['Finished nominal Scenario Calculation. Computation time: ', num2str(round(nomScenTime / 3600, 2)), 'h \n'],param,'info');

refVol = [2 5 50 95 98];
refGy = linspace(0,max(resultGUInomScen.(pln.bioParam.quantityVis)(:)),6);

resultGUInomScen.dvh = matRad_calcDVH(cst,resultGUInomScen.(pln.bioParam.quantityVis),'cum');
dvhPoints            = resultGUInomScen.dvh(1).doseGrid;
nomQi                = matRad_calcQualityIndicators(cst,pln,resultGUInomScen.(pln.bioParam.quantityVis),refGy,refVol,param);

resultGUInomScen.qi  = nomQi;
resultGUInomScen.cst = cst;

nomDose                = resultGUInomScen.(pln.bioParam.quantityVis)(param.subIx);
nomDoseLin             = single(reshape(nomDose,[],1));

nomScen = MatRadScenario(nomDoseLin, param.subIx, pln.bioParam.quantityVis, 1, [0 0 0], 0, 0, size(resultGUInomScen.(pln.bioParam.quantityVis)), NaN);

% only show warnings and disable waitbars and figures
param.logLevel = 3;

%% perform parallel sampling
scenContainer = cell(pln.multScen.totNumScen,1);
if FlagParallToolBoxLicensed
   % Create parallel pool on cluster
   p = gcp(); % If no pool, create new one.
   
   if isempty(p)
       poolSize = 1;
   else
       poolSize = p.NumWorkers;
   end
   % rough estimate of total computation time
   totCompTime = ceil(size(pln.multScen.scenForProb,1) / poolSize) * nomScenTime * 1.35;
   fprintf(['Approximate Total calculation time: ', num2str(round(totCompTime / 3600, 2)), ...
                      'h. Estimated finish: ', datestr(datetime('now') + seconds(totCompTime)), '\n']);
   
   if exist('parfor_progress', 'file') == 2
      FlagParforProgressDisp = true;
      parfor_progress(pln.multScen.totNumScen);  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
   else
      fprintf('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n');
      FlagParforProgressDisp = false;
   end
   
   parfor i = 1:pln.multScen.totNumScen
          
          % create nominal scenario
          plnSamp          = pln;
          plnSamp.multScen = matRad_multScen(ct,'nomScen'); 
          % pick the i-th scenario and save into plnSamp
          plnSamp.multScen.scenForProb         = pln.multScen.scenForProb(i,:);
          plnSamp.multScen.relRangeShift       = pln.multScen.scenForProb(i,5);
          plnSamp.multScen.absRangeShift       = pln.multScen.scenForProb(i,4);
          plnSamp.multScen.isoShift            = pln.multScen.scenForProb(i,1:3);
          plnSamp.multScen.totNumShiftScen     = 1;
          plnSamp.multScen.totNumRangeScen     = 1;
          plnSamp.multScen.numOfCtScen         = 1;
          plnSamp.multScen.scenMask            = 1;
          plnSamp.multScen.linearMask          = 1;
          plnSamp.multScen.scenProb            = 1;
          
          resultSamp                 = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose                = resultSamp.(pln.bioParam.quantityVis)(param.subIx);
          sampledDoseLin             = single(reshape(sampledDose,[],1));
          scenContainer{i}           = MatRadScenario(sampledDoseLin, param.subIx, pln.bioParam.quantityVis, 1, plnSamp.multScen.isoShift, ...
                plnSamp.multScen.relRangeShift, plnSamp.multScen.absRangeShift, size(resultSamp.(pln.bioParam.quantityVis)), 1);
          
          % caSampRes(i).dvh = matRad_calcDVH(cst,resultSamp.(pln.bioParam.quantityOpt),'cum',dvhPoints);
          % caSampRes(i).qi  = matRad_calcQualityIndicators(cst,pln,resultSamp.(pln.bioParam.quantityOpt),refGy,refVol,param);
          
          if FlagParforProgressDisp
            parfor_progress;
          end
   end

   if FlagParforProgressDisp
      parfor_progress(0);
   end

else
%% perform seriel sampling
% rough estimate of total computation time
totCompTime = size(pln.multScen.scenForProb,1) * nomScenTime * 1.1;
fprintf(['Approximate Total calculation time: ', num2str(round(totCompTime / 3600, 2)), ...
                        'h. Estimated finish: ', datestr(datetime('now') + seconds(totCompTime)), '\n']);
    
    for i = 1:pln.multScen.totNumScen
       
          % create nominal scenario
          plnSamp          = pln;
          plnSamp.multScen = matRad_multScen(ct,'nomScen'); 
          % pick the i-th scenario and save into plnSamp
          plnSamp.multScen.scenForProb         = pln.multScen.scenForProb(i,:);
          plnSamp.multScen.relRangeShift       = pln.multScen.scenForProb(i,5);
          plnSamp.multScen.absRangeShift       = pln.multScen.scenForProb(i,4);
          plnSamp.multScen.isoShift            = pln.multScen.scenForProb(i,1:3);
          plnSamp.multScen.totNumShiftScen     = 1;
          plnSamp.multScen.totNumRangeScen     = 1;
          plnSamp.multScen.numOfCtScen         = 1;
          plnSamp.multScen.scenMask            = 1;
          plnSamp.multScen.linearMask          = 1;
          plnSamp.multScen.scenProb            = 1;
          
          resultSamp                 = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose                = resultSamp.(pln.bioParam.quantityVis)(param.subIx);
          sampledDoseLin             = single(reshape(sampledDose,[],1));
          scenContainer{i}           = MatRadScenario(sampledDoseLin, param.subIx, pln.bioParam.quantityVis, 1, plnSamp.multScen.isoShift, ...
                plnSamp.multScen.relRangeShift, plnSamp.multScen.absRangeShift, size(resultSamp.(pln.bioParam.quantityVis)), 1);
          
          % caSampRes(i).dvh = matRad_calcDVH(cst,resultSamp.(pln.bioParam.quantityOpt),'cum',dvhPoints);
          % caSampRes(i).qi  = matRad_calcQualityIndicators(cst,pln,resultSamp.(pln.bioParam.quantityOpt),refGy,refVol,param);

          matRad_progress(i, pln.multScen.totNumScen)
    end
    
end

%% add subindices
pln.subIx        = param.subIx;

%% account for fractionation
treatmentSimulation = MatRadSimulation(pln.bioParam.quantityVis, nomScen, ct, cst, pln, pln.numOfFractions);
for i = 1:numel(scenContainer)
    treatmentSimulation.initNewScen(scenContainer{i})
end
end
