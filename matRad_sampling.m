function [mRealizations,stats, pln, resultCubes]  = matRad_sampling(ct,stf,cst,pln,w,structSel, param)
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
pln.numOfSamples  = 4;


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

stats       = cell(pln.numOfSamples,1);

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
          mRealizations(:,i)     = single(reshape(sampledDose,[],1));
          resultQI              = matRad_calcQualityIndicators(resultSamp,cst,plnSamp,param);
          stats{i,1}            = resultQI.QI;
          
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
          resultQI              = matRad_calcQualityIndicators(resultSamp,cst,plnSamp,param);
          stats{i,1}            = resultQI.QI;
          
          waitbar(i/pln.numOfSamples);

    end
    
    close(h)
end

%% add nominal scenario
resultGUInominal          = matRad_calcDoseDirect(ct,stf,pln,cst,w,param);
resultCubes.resultNominal = resultGUInominal.(pln.bioParam.quantityOpt);      
        
%% add subindices
pln.multScen.subIx        = param.subIx;
resultCubes.subIx         = param.subIx;

end
