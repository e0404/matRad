function dij = matRad_cleanDijScenarios(dij,pln,cst)
% matRad function to clean up dij in scenarios after dose calculation to 
%   remove dose influence for voxels outside of segmentations for every ct
%   scenario
%
% call
%   dij = matRad_calcPhotonDose(dij,pln,cst)
%
% input
%   dij:            matRad dij struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%
% output
%   dij:            "clean" matRad dij struct
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:pln.multScen.numOfCtScen
   
   % generate index set to erase
   tmpIx = [];
   for j = 1:size(cst,1)
      tmpIx = unique([tmpIx; cst{j,4}{i}]);
   end
   ix = setdiff(1:dij.doseGrid.numOfVoxels,tmpIx);
   
   for j = 1:pln.multScen.totNumShiftScen
      for k = 1:pln.multScen.totNumRangeScen
         
         if pln.multScen.scenMask(i,j,k)
            
            dij.physicalDose{i,j,k}(ix,:)      = 0;
            
            if isfield(dij,'mLETDose')
               dij.mLETDose{i,j,k}(ix,:)      = 0;
            end
            
            if pln.bioParam.bioOpt
               dij.mAlphaDose{i,j,k}(ix,:)    = 0;
               dij.mSqrtBetaDose{i,j,k}(ix,:) = 0;
            end
            
         end
         
      end
   end
end

end

