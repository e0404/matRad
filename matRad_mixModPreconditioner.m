function dij = matRad_mixModPreconditioner(dij) 
% Dose influence preconditioner for mixed modality plans
%
%
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


% check if its is a mixed mod dij exists 
if ~isfield(dij,'original_Dijs')
    return;
end

% check to identify what is the preconditioning value 

for i = 1: numel(dij.original_Dijs)
    preconW(i) = max(mean(dij.original_Dijs{i}.physicalDose{1},1));
end
preconW = round(preconW./min(preconW));
dij.preconW = 1./preconW;

% make change to all dij type structures in one 
for mod = 1 : numel(dij.original_Dijs)
    fieldNames = fieldnames(dij.original_Dijs{mod});

    for i = 1 : numel(fieldNames)
        if iscell(dij.original_Dijs{mod}.(fieldNames{i}))
            
            for j = 1 : numel(dij.original_Dijs{mod}.(fieldNames{i})) 
                dij.original_Dijs{mod}.(fieldNames{i}){j} = bsxfun(@times,dij.original_Dijs{mod}.(fieldNames{i}){j}, dij.preconW(mod)); 
            end
        end
    end
end

end