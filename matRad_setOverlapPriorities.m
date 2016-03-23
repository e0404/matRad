function cst = matRad_setOverlapPriorities(cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to handle overlap priorities during fluence optimizaiton and
% dose calculation. If you have overlapping volumes of interest you need to
% inform matrad to which volume(s) the intersection voxels belong
% 
% call
%   cst = matRad_considerOverlap(cst)
%
% input
%   cst:    cst file
%
% output
%   cst:    updated cst file considering overlap priorities
%
% References
%   -
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


    % consider VOI priorities
    for i = 1:size(cst,1)
         
        idx = cst{i,4};          
        
        for k = 1:size(cst,1)
            if cst{k,5}.Priority < cst{i,5}.Priority && ~(i==k)
                % remove indices from VOI with higher priority from current VOI
                idx = setdiff(idx,cst{k,4});
            end
        end
        
        cst{i,4} = idx;
        
        if isempty(cst{i,4}) && ~isempty(cst{i,6})
            warning([cst{i,2} ': Objective(s) for inverse planning defined ' ...
                 'but structure overlapped by structure with higher overlap priority.' ...
                 'Objective(s) will not be considered during optimization']); 
        end
         
    end

end

