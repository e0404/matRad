function [sampRes, sampDose, pln, nominalScenario]  = matRad_sampling(ct,stf,cst,pln,w,structSel, multScen, param)
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
%   sampDose:  matrix depticting the sampling results in V x pln.numSamples
%   stats          cell array denoting dose statistics
%   meanCube       expected dose distribution
%   stdCube        standard deviation cube
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

pln.sampling      = true;
pln.robOpt        = false;

% calculate RBExDose
if strcmp(pln.radiationMode, 'protons')
    pln.bioOptimization = 'constRBE_RBExD';
else strcmp(pln.radiationMode, 'carbon')
    pln.bioOptimization = 'LEM_RBExD';
end
pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);

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


% save nonSampling pln for nominal scenario calculation
plnNominal = pln;
plnNominal.multScen.numOfRangeShiftScen = 1;
if exist('multScen','var')
    pln = matRad_setPlanUncertainties(ct,pln, multScen, param);
else
    pln = matRad_setPlanUncertainties(ct,pln, [], param);
end


matRad_dispToConsole(['Using ' num2str(pln.multScen.numOfScen) 'samples in total \n'],param,'info')

V = [];
% define voxels for sampling
if ~exist('structSel', 'var') || sum(size(structSel)) == 0
    V = [cst{:,4}];
else
    for i=1:size(cst,1)
        for j = 1:numel(structSel)
            if strcmp(structSel{j}, cst{i,2})
                V = [V cst{i, 4}];
            end
        end
    end
end

param.subIx = unique(vertcat(V{:}));

% disable structures which are not completely in subIx
for i = 1:size(cst,1)
    if ~all(ismember(cst{i,4}{1}, param.subIx))
        cst{i,5}.Visible = false;
    end
end

% define variable for storing scenario doses
sampDose = single(zeros(numel(param.subIx),pln.numOfSamples,1));
StorageInfo  = whos('sampDose');
matRad_dispToConsole(['matRad: Realizations variable will need: ' num2str(StorageInfo.bytes/1e9) ' GB \n'],param,'info');


% only show warnings
param.logLevel = 3;

% check if parallel toolbox is installed and license can be checked out
try
   ver('distcomp')                 
   FlagParallToolBoxLicensed  = license('test','Distrib_Computing_Toolbox'); 
catch
   FlagParallToolBoxLicensed  = false;
end

%% calculate nominal scenario
nomScenTimer = tic;
nominalScenario          = matRad_calcDoseDirect(ct,stf,plnNominal,cst,w,param);
nomScenTime = toc(nomScenTimer);
fprintf(['Finished nominal Scenario Calculation. Computation time: ', num2str(round(nomScenTime / 3600, 2)), 'h \n']);

refVol = [2 5 50 95 98];
refGy = linspace(0,max(nominalScenario.(pln.bioParam.quantityOpt)(:)),6);

nominalScenario.dvh = matRad_calcDVH(cst,nominalScenario.(pln.bioParam.quantityOpt),'cum');
dvhPoints = nominalScenario.dvh(1).doseGrid;
nomQi = matRad_calcQualityIndicators(cst,pln,nominalScenario.(pln.bioParam.quantityOpt),refGy,refVol,param);

nominalScenario.qi = nomQi;
nominalScenario.cst = cst;
%% perform parallel sampling

if FlagParallToolBoxLicensed
   % Create parallel pool on cluster
   p = gcp('nocreate'); % If no pool, do not create new one.
   
   % rough estimate of total computation time
   totCompTime = size(pln.multScen.scenForProb,1) * nomScenTime * 1.25 / p.NumWorkers;
   fprintf(['Approximate Total calculation time: ', num2str(round(totCompTime / 3600, 2)), ...
                      'h. Estimated finish: ', datestr(datetime('now') + seconds(totCompTime)), '\n']);
   
   if exist('parfor_progress', 'file') == 2
      FlagParforProgressDisp = true;
      parfor_progress(pln.numOfSamples);  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
   else
      Console('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n',param,'warning');
      FlagParforProgressDisp = false;
   end
  
   parfor i = 1:size(pln.multScen.scenForProb,1)
          
          plnSamp               = pln;
          % pick the i-th scenario and save into plnSamp
          plnSamp.multScen.scenForProb         = pln.multScen.scenForProb(i,:);
          plnSamp.multScen.relRangeShift       = pln.multScen.scenForProb(i,5);
          plnSamp.multScen.absRangeShift       = pln.multScen.scenForProb(i,4);
          plnSamp.multScen.isoShift            = pln.multScen.scenForProb(i,1:3);
          plnSamp.multScen.numOfShiftScen      = 1;
          plnSamp.multScen.numOfRangeShiftScen = 1;
          plnSamp.multScen.numOfCtScen         = 1;
          plnSamp.multScen.scenMask            = 1;
          plnSamp.multScen.linearMask          = 1;
          plnSamp.multScen.scenProb            = 1;
          
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          sampDose(:,i)         = single(reshape(sampledDose,[],1));
          sampRes(i).bioParam      = pln.bioParam;
          sampRes(i).relRangeShift = plnSamp.multScen.relRangeShift;
          sampRes(i).absRangeShift = plnSamp.multScen.absRangeShift;
          sampRes(i).isoShift      = plnSamp.multScen.isoShift;
          
          sampRes(i).dvh = matRad_calcDVH(cst,resultSamp.(pln.bioParam.quantityOpt),'cum',dvhPoints);
          sampRes(i).qi = matRad_calcQualityIndicators(cst,pln,resultSamp.(pln.bioParam.quantityOpt),refGy,refVol,param);
          
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
    
    for i = 1:pln.numOfSamples
       
          plnSamp               = pln;
          % pick the i-th scenario and save into plnSamp
          plnSamp.multScen.scenForProb         = pln.multScen.scenForProb(i,:);
          plnSamp.multScen.relRangeShift       = pln.multScen.scenForProb(i,5);
          plnSamp.multScen.absRangeShift       = pln.multScen.scenForProb(i,4);
          plnSamp.multScen.isoShift            = pln.multScen.scenForProb(i,1:3);
          plnSamp.multScen.numOfShiftScen      = 1;
          plnSamp.multScen.numOfRangeShiftScen = 1;
          plnSamp.multScen.numOfCtScen         = 1;
          plnSamp.multScen.scenMask            = 1;
          plnSamp.multScen.linearMask          = 1;
          plnSamp.multScen.scenProb            = 1;
          
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          sampDose(:,i)         = single(reshape(sampledDose,[],1));
          sampRes(i).bioParam      = pln.bioParam;
          sampRes(i).relRangeShift = plnSamp.multScen.relRangeShift;
          sampRes(i).absRangeShift = plnSamp.multScen.absRangeShift;
          sampRes(i).isoShift      = plnSamp.multScen.isoShift;
          
          sampRes(i).dvh = matRad_calcDVH(cst,resultSamp.(pln.bioParam.quantityOpt),'cum',dvhPoints);
          sampRes(i).qi = matRad_calcQualityIndicators(cst,pln,resultSamp.(pln.bioParam.quantityOpt),refGy,refVol,param);
          matRad_progress(i, pln.numOfSamples)
    end
    

end
%% add subindices
pln.multScen.subIx        = param.subIx;

end
