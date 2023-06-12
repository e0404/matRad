function stf = matRad_stfWrapper(ct,cst,plnJO)
% wrapper to handle structure field consistency and to cmobine multimodal
% stf structures in to one stf centipede
% 
% call
%   stf = matRad_stfWrapper(ct,cst,plnJO)
%
% input
%   ct:         matRad ct stuct 
%   cst:        matrad cst struct
%   plnJO:      Joint opt pln struct 
%
% output
%   stf:        matRad stf stuct merged for all modalities 
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   stf = [];
   fields = [];   
   
   %Generate stf for every modality
   for k=1:length(plnJO.propStf)

      currStf = matRad_generateStf(ct,cst,plnJO.originalPlans(k));
      fields = [fields; fieldnames(currStf)];
      stf = [stf,{currStf}];
      
   end
   
   %Add fields to the single stf structures when missing in order to make
   %them compatible
   totalFields = unique(fields);
   for k=1:length(stf)
      isPlanField = find(~isfield(stf{k},totalFields));
      if any(isPlanField)
         for m=[isPlanField]
            stf{k} = setfield(stf{k},{1},totalFields{m},[]);
         end
      end
   end
   stf = [stf{:}];
end
