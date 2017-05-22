function [Realizations,Stats,meanCube,stdCube]  = matRad_randomSampling(ct,stf,cst,pln,w)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_randomSampling enables sampling multiple treatment scenarios
% 
% call
%   [cst,pln] = matRad_setPlanUncertainties(ct,cst,pln)
%
% input
%   ct:             ct cube
%   pln:            matRad plan meta information struct
%
% output
%   pln:            matRad's plan meta information struct including a sub-structure 
%                   pln.multScen holding information about multiple treatment plan scenarios
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

global LogLevel;

if ~isfield(pln,'sampling')
   error('matRad: pln.samling subfield does not exist');
else
   if ~isfield(pln,'numSampling')
      pln.numSampling  = 30; % default number of samples
   end
end

meanCube    = zeros(ct.cubeDim);
stdCube     = zeros(ct.cubeDim);
Stats       = 0;

% define voxels for sampling
V = [cst{:,4}];
param.subIx = unique(vertcat(V{:}));

% define variable for storing scenario doses
Realizations = single(zeros(numel(param.subIx),pln.numSampling,1));
StorageInfo  = whos('Realizations');
fprintf(['matRad: Realizations variable will need: ' num2str(StorageInfo.bytes/1e9) ' GB \n']);


% if sampling is disabled stop this function call
if pln.sampling == false
    return;
end

% only show warning an
LogLevel = 3;

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
      warning('matRad: Consider downloading parfor_progress function from the matlab central fileexchange to get feedback from loop.');
      FlagParforProgressDisp = false;
   end
  
   plnTot               = matRad_setPlanUncertainties(ct,pln);
   
   parfor i = 1:pln.numSampling
    
          % ToDo pick the i-th scenario and save into plnSamp
          plnSamp               = plnTot;
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          SampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          Realizations(:,i)     = single(reshape(SampledDose,[],1));
          resultQI              = matRad_calcQualityIndicators(resultSamp,cstSampling,plnSampling,param);
          Stats{i}              = resultQI.QI;
          
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
    Stats = cell(pln.numSampling,1);
    
    plnTot               = matRad_setPlanUncertainties(ct,pln);
    
    for i = 1:pln.numSampling
       
           % ToDo pick the i-th scenario and save into plnSamp
          plnSamp               = plnTot; 
          resultSamp            = matRad_calcDoseDirect(ct,stf,plnSamp,cst,w,param);
          SampledDose           = resultSamp.(pln.bioParam.quantityOpt)(param.subIx);
          Realizations(:,i)     = single(reshape(SampledDose,[],1));
          resultQI              = matRad_calcQualityIndicators(resultSamp,cstSampling,plnSampling,param);
          Stats{i}              = resultQI.QI;
          
          waitbar(i/pln.numSampling);

    end
    
    close(h)
end

meanCube    = zeros(ct.cubeDim);
stdCube     = zeros(ct.cubeDim);
meanCube(V) = mean(Realizations,2);
stdCube(V)  = std(Realizations,1,2);


end

