function  [dij] = matRad_estimateParticleMean(cst,pln,dij)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on the calculated dose scenarios in the dij structure,
% the expectation value will be estimated and stored in the dij struct
% variance in the cst.
% 
% call
%    [dij] = matRad_calcProbParticleDoseEstimate(cst,pln,dij)
%
% input
%
%   cst:            matRad cst struct
%   pln:            matRad plan meta information struct
%   dij:            matRad dij struct
%
% output
%   dij:            matRad dij struct
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
% 

if ~pln.bioParam.bioOpt
   fNames = {'physicalDose'};
else
   fNames = {'physicalDose','mAlphaDose','mSqrtBetaDose'};
end

ixDij = find(~cellfun(@isempty, dij.physicalDose))';

for i = 1:numel(fNames)
   % create expected ij structure
   dij.([fNames{1,i} 'Exp']){1} = spalloc(prod(dij.dimensions),dij.totalNumOfBixels,1);
   % add up sparse matrices - should possess almost same sparsity pattern
   for j = 1:pln.multScen.totNumScen
      dij.([fNames{1,i} 'Exp']){1} = dij.([fNames{1,i} 'Exp']){1} + dij.([fNames{1,i}]){ixDij(j)} .* pln.multScen.scenProb(j);
   end
    
end



end