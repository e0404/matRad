function  [cst,resultGUI] = matRad_estimateParticleVar(cst,pln,dij,resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the integral variance based on the calculated dose scenarios in the dij structure,
% 
% call
%    [cst] = matRad_calcProbParticleDoseEstimate(cst,pln,dij)
%
% input
%
%   cst:            matRad cst struct
%   pln:            matRad plan meta information struct
%   dij:            matRad dij struct
%
% output
%   cst:            matRad cst struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
% 
warning('code is not thoroughly tested');

if ~pln.bioParam.bioOpt
   fNames = {'physicalDose'};
else
   fNames = {'mAlphaDose'};
   matRad_dispToConsole('Only variance in mAlphaDose is considered',[],'warning');
end

%% calculate the standard deviation of the dose cube

voiIx = find(~cellfun(@isempty, cst(:,6)))'; 

for i = 1:numel(fNames)
   for j = voiIx   
      cst{j,6}(1).mOmega = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
      
      % loop over scenarios and calculate the integral variance of each
      % spot combination; bio bio optimization only consider std in the
      % linear part of the biological effect
      for k = 1:pln.multScen.totNumScen
         cst{j,6}(1).mOmega = cst{j,6}(1).mOmega + ...
            ((dij.(fNames{i,1}){ixDij(k)}(cst{j,4}{1},:)' * dij.(fNames{i,1}){ixDij(k)}(cst{j,4}{1},:)) * pln.multScen.scenProb(k));
      end
      cst{j,6}(1).mOmega = cst{j,6}(1).mOmega - (dij.([fNames{i,1} 'Exp']){1}(cst{j,4}{1},:)' * dij.([fNames{i,1} 'Exp']){1}(cst{j,4}{1},:));
   end
end

%% To do: calculate variance and add it to resultGUI


end