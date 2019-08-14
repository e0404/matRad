function cst = matRad_computeVoiContoursWrapper(cst,ct)
% matRad computation of VOI contours if not precomputed
% 
% call
%   cst = matRad_computeVoiContoursWrapper(ct,cst)
% 
% input:
%   cst:        matRad cst struct
%   ct:         matRad ct struct
% 
% output:
%   cst:        matRad cst struct with VOI contours
% 
% References
%   -
% 
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

if size(cst,2) < 7
    cst = matRad_computeVoiContours(ct,cst);
else
    for i = 1:size(cst,1)
        if isempty(cst{i,7})
            cst = matRad_computeVoiContours(ct,cst);
            break
        end
    end
end
