function resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst,w,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dose calculation wrapper bypassing dij calculation
%
% call
%   resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   w:          optional (if no weights available in stf): bixel weight
%               vector
%
% output
%   resultGUI:  matRad result struct
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
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
else
   param.subIx          = [];
   param.logLevel       = 1;
end

param.calcDoseDirect = true;

% check if weight vector is available, either in function call or in stf - otherwise dose calculation not possible
if ~exist('w','var') && ~isfield([stf.ray],'weight')
   error('No weight vector available. Please provide w or add info to stf')
end

% copy bixel weight vector into stf struct
if exist('w','var')
   if sum([stf.totalNumOfBixels]) ~= numel(w)
      error('weighting does not match steering information')
   end
   counter = 0;
   for i = 1:size(stf,2)
      for j = 1:stf(i).numOfRays
         for k = 1:stf(i).numOfBixelsPerRay(j)
            counter = counter + 1;
            stf(i).ray(j).weight(k) = w(counter);
         end
      end
   end
else % weights need to be in stf!
   w = NaN*ones(sum([stf.totalNumOfBixels]),1);
   counter = 0;
   for i = 1:size(stf,2)
      for j = 1:stf(i).numOfRays
         for k = 1:stf(i).numOfBixelsPerRay(j)
            counter = counter + 1;
            w(counter) = stf(i).ray(j).weight(k);
         end
      end
   end
end

% dose calculation
[cst,dij] = matRad_calcDose(ct,stf,pln,cst,param);

% calc resulting dose
if pln.multScen.totNumScen == 1
   % calculate cubes; use uniform weights here, weighting with actual fluence
   % already performed in dij construction
   resultGUI    = matRad_calcCubes(ones(pln.propStf.numOfBeams,1),dij,cst);
   
   % calc individual scenarios
else
   
   Cnt          = 1;
   ixForOpt     = find(~cellfun(@isempty, dij.physicalDose))';
   for i = ixForOpt
      tmpResultGUI = matRad_calcCubes(ones(pln.propStf.numOfBeams,1),dij,cst,i);
      resultGUI.([pln.bioParam.quantityVis '_' num2str(Cnt,'%d')]) = tmpResultGUI.(pln.bioParam.quantityVis);
      Cnt = Cnt + 1;
   end
end

% remember original fluence weights
resultGUI.w  = w;




