function [mRealizations,stats, cst, pln, resultCubes,nominalScenario]  = matRad_sampling(ct,stf,cst,pln,w,structSel, param)
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
%   mRealizations:  matrix depticting the sampling results in V x pln.numSamples
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
pln.numOfSamples  = 108;


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


if ~isfield(pln,'numOfSamples')
   pln.numOfSamples  = 20; % default number of samples
end
matRad_dispToConsole(['Using samples: ' num2str(pln.numOfSamples) ' in total \n'],param,'info')

stats       = cell(pln.numOfSamples,2);
% since parfor does not allow different calling
dvhTemp     = cell(pln.numOfSamples,1);

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

% define variable for storing scenario doses
mRealizations = single(zeros(numel(param.subIx),pln.numOfSamples,1));
StorageInfo  = whos('mRealizations');
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
nominalScenario          = matRad_calcDoseDirect(ct,stf,pln,cst,w,param);

nominalScenario.cst = cst;
nominalScenario.dvh = matRad_calcDVH(cst,nominalScenario.(pln.bioParam.quantityOpt),'cum');

refVol = [2 5 98 95];
refGy = linspace(0,max(nominalScenario.(pln.bioParam.quantityOpt)(:)),6);
nomQi = matRad_calcQualityIndicators(cst,pln,nominalScenario.(pln.bioParam.quantityOpt),refGy,refVol);
nominalScenario.qi = nomQi;
for i = 1:size(nominalScenario.cst,1)
    nominalScenario.cst{i,8} = cell(1,1);
    nominalScenario.cst{i,9} = cell(1,1);
    
    nominalScenario.cst{i,8}{1} = nominalScenario.dvh{i};
    nominalScenario.cst{i,9}{1} = nomQi{i};
end  
dvhPoints = nominalScenario.dvh{1}(1,:);

%% perform parallel sampling

if FlagParallToolBoxLicensed
   % Create parallel pool on cluster
   p = gcp('nocreate'); % If no pool, do not create new one.
   
   if exist('parfor_progress', 'file') == 2
      FlagParforProgressDisp = true;
      parfor_progress(pln.numOfSamples);  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
   else
      matRad_dispToConsole('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.\n',param,'warning');
      FlagParforProgressDisp = false;
   end
  
   pln               = matRad_setPlanUncertainties(ct,pln);
   
   parfor i = 1:pln.numOfSamples
          
          plnSamp               = pln;
          % pick the i-th scenario and save into plnSamp
          plnSamp.multScen.scenForProb     = pln.multScen.scenForProb(i,:);
          plnSamp.multScen.relRangeShift   = pln.multScen.relRangeShift(i);
          plnSamp.multScen.absRangeShift   = pln.multScen.absRangeShift(i);
          plnSamp.multScen.numOfShiftScen  = 1;
          plnSamp.multScen.numOfRangeShift = 1;
          plnSamp.multScen.numOfCtScen     = 1;
          plnSamp.multScen.scenMask        = 1;
          plnSamp.multScen.linearMask      = 1;
          plnSamp.multScen.scenProb        = 1;
          
          
          % forward dose calculation
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          mRealizations(:,i)    = single(reshape(sampledDose,[],1));
          
          dvh{i} = matRad_calcDVH(cst,resultSamp.(pln.bioParam.quantityOpt),'cum',dvhPoints);
          qi{i} = matRad_calcQualityIndicators(cst,pln,resultSamp.(pln.bioParam.quantityOpt),refGy,refVol);
          
          % resultQI              = matRad_calcQualityIndicators(resultSamp,cst,plnSamp,param);
          % stats{i,1}            = resultQI.QI;
          % [~,dvhTemp{i,1}]      = calcDVHdirect(cst, resultSamp.(pln.bioParam.quantityOpt), 1, resultCubes.dvhPoints);
          
          if FlagParforProgressDisp
            parfor_progress;
          end
   end

   if FlagParforProgressDisp
      parfor_progress(0);
   end

else
%% perform seriel sampling   
    h = waitbar(0,'Sampling Scenario ...');
    stats = cell(pln.numOfSamples,1);
    
    pln               = matRad_setPlanUncertainties(ct,pln);
    
    for i = 1:pln.numOfSamples
       
          plnSamp               = pln;
          % pick the i-th scenario and save into plnSamp
          plnSamp.multScen.scenForProb     = pln.multScen.scenForProb(i,:);
          plnSamp.multScen.relRangeShift   = pln.multScen.relRangeShift(i);
          plnSamp.multScen.absRangeShift   = pln.multScen.absRangeShift(i);
          plnSamp.multScen.numOfShiftScen  = 1;
          plnSamp.multScen.numOfRangeShift = 1;
          plnSamp.multScen.numOfCtScen     = 1;
          plnSamp.multScen.scenMask        = 1;
          plnSamp.multScen.linearMask      = 1;
          plnSamp.multScen.scenProb        = 1;
          
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          mRealizations(:,i)     = single(reshape(sampledDose,[],1));
          
          dvh{i} = matRad_calcDVH(cst,sampledDose,'cum',dvhPoints);
          qi{i} = matRad_calcQualityIndicators(cst,pln,sampleDose,refGy,refVol);
          
          waitbar(i/pln.numOfSamples);

    end
    
    close(h)
end

% reassing dvh to stats structure
for i = 1:size(nominalScenario.cst,1)
    cst{i,8} = cell(pln.numOfSamples,1);
    cst{i,9} = cell(pln.numOfSamples,1);
    for j = 1:pln.numOfSamples
        cst{i,8}{j} = dvh{j}{i};
        cst{i,9}{j} = qi{j}{i};
    end  
end
        
%% add subindices
pln.multScen.subIx        = param.subIx;
resultCubes.subIx         = param.subIx;

end
