function [realizations,stats,meanCube,stdCube]  = matRad_randomSampling(ct,stf,cst,pln,w,param)
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
%   realizations:  matrix depticting the sampling results in V x pln.numSamples
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


if ~isfield(pln,'numSampling')
   pln.numSampling  = 20; % default number of samples
   matRad_dispToConsole(['Using default number of samples: ' num2str(pln.numSampling) '\n'],param,'info')
end

meanCube    = zeros(ct.cubeDim);
stdCube     = zeros(ct.cubeDim);
stats       = cell(pln.numSampling,1);

% define voxels for sampling
V = [cst{:,4}];
param.subIx = unique(vertcat(V{:}));

% define variable for storing scenario doses
realizations = single(zeros(numel(param.subIx),pln.numSampling,1));
StorageInfo  = whos('realizations');
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
      parfor_progress(pln.numSampling);  % http://de.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor
   else
      matRad_dispToConsole('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from parfor loop.',param,'warning');
      FlagParforProgressDisp = false;
   end
  
   plnTot               = matRad_setPlanUncertainties(ct,pln);
   
   parfor i = 1:pln.numSampling
    
          % ToDo pick the i-th scenario and save into plnSamp
          plnSamp               = plnTot;
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          realizations(:,i)     = single(reshape(sampledDose,[],1));
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
    stats = cell(pln.numSampling,1);
    
    plnTot               = matRad_setPlanUncertainties(ct,pln);
    
    for i = 1:pln.numSampling
       
           % ToDo pick the i-th scenario and save into plnSamp
          plnSamp               = plnTot; 
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          sampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          realizations(:,i)     = single(reshape(sampledDose,[],1));
          resultQI              = matRad_calcQualityIndicators(resultSamp,cst,plnSamp,param);
          stats{i,1}            = resultQI.QI;
          
          waitbar(i/pln.numSampling);

    end
    
    close(h)
end

meanCube              = zeros(ct.cubeDim);
stdCube               = zeros(ct.cubeDim);
meanCube(param.subIx) = mean(realizations,2);
stdCube(param.subIx)  = std(realizations,1,2);


end

