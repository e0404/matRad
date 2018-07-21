function [cst,dijExp] = matRad_calcDose(ct,stf,pln,cst,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dose influence calculation wrapper
%
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   param:          (optional) structure defining additional parameter
%                   param.calcDoseDirect boolian switch to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
%
% output
%   dij:            matRad dij struct
%   cst:            matRad cst struct
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('APM')

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

%% set covariance matrices
pln.multScen.setCovarianceMatrix(ct,cst,pln,stf);

%% set nominal dose calculation to true
param.nomDoseCalc = true;    % determines if the nominal or the expected dose should be calcuated

% calculate scenario doses or nominal doses
if strcmp(pln.radiationMode,'photons')
   dij = matRad_calcPhotonDose(ct,stf,pln,cst,param);
   %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
   dij = matRad_calcParticleDose(ct,stf,pln,cst,param);
end

%% start calculating the expected dose influence matrices
matRad_dispToConsole('Calculating probabilistic quantities for optimization ...\n',param,'info');

% get fieldnames
if ~pln.bioParam.bioOpt
   fNames = {'physicalDose'};
else
   fNames = {'physicalDose','mAlphaDose','mSqrtBetaDose'};
end

if  isfield(pln,'propDoseCalc')
    if isfield(pln.propDoseCalc,'calcLET')
        if pln.propDoseCalc.calcLET
            fNames{1,end+1} = 'mLETDose';
        end
    end
end

% perform either analytical calculation, estimation or do nothing
switch pln.multScen.TYPE
   
   % perform analytical first moment calculation
   case 'apmScen'
      
      param.nomDoseCalc = false;
      options           = matRad_probOptions(ct,cst,pln);
      pln.probOpt       = options.probOpt;
      
      % calcualted expected dose influence matrices
      dijExp = matRad_calcParticleDose(ct,stf,pln,cst,param);
      
      suffix = 'Exp';
      % organize fieldnames
      for i = 1:numel(fNames)
         % create new Exp field and set it
         [dijExp.([fNames{1,i} suffix])] = dijExp.([fNames{1,i}]);
         % overide existing field with nominal dose influence data field
         dijExp.([fNames{1,i}]) = dij.([fNames{1,i}]);
      end
      
      % find VOI indicies with objective or constraint and create a field in the cst for the integral variance
      cst = matRad_createOmegaPlaceHolder(cst,dij.totalNumOfBixels);
      
      % estimate first moment based on discrete samples
   case  {'wcScen','impScen'}
      
      [dijExp] = matRad_estimateParticleMean(cst,pln,dij);
      cst      = matRad_createOmegaPlaceHolder(cst,dijExp.totalNumOfBixels);
      
      % otherwise create empty placeholders for expected ij matrices
   otherwise
      
      dijExp = dij;
      clear 'dij';
      
      if param.calcDoseDirect
         numCol =  dijExp.numOfBeams;
      else
         numCol = dijExp.totalNumOfBixels;
      end
      
      for i = 1:numel(fNames)
         dijExp.([fNames{1,i} 'Exp']){1} = spalloc(prod(dijExp.dimensions),numCol,1);
      end
      
      % find VOI indicies with objective or constraint and create a field in the cst for the integral variance
      cst = matRad_createOmegaPlaceHolder(cst,dijExp.totalNumOfBixels);
      
end %eof switch



end % eof function


function cst = matRad_createOmegaPlaceHolder(cst,numBixels)
voiIx = find(~cellfun(@isempty, cst(:,6)))';
for i = voiIx
   cst{i,6}(1).mOmega = spalloc(numBixels,numBixels,1);
end
end


