function cst = matRad_computeVoiContours(ct,cst)
% matRad function that computes all VOI contours
%
% call
%   cst = matRad_computeVoiContours(ct,cst)
%
% input
%   ct  matRad ct struct
%   cst matRad cst struct
%
% output
%   cst the new cst with the column containing the precomputed contours
%
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

mask = zeros(ct.cubeDim); % create zero cube with same dimeonsions like dose cube
for s = 1:size(cst,1)
    cst{s,7} = cell(max(ct.cubeDim(:)),3);
    mask(:) = 0;
    mask(cst{s,4}{1}) = 1;    
    for slice = 1:ct.cubeDim(1)
        if any(any(mask(slice,:,:) > 0))
             cst{s,7}{slice,1} = contourc(squeeze(mask(slice,:,:)),.5*[1 1]);
        end
    end
    for slice = 1:ct.cubeDim(2)
        if any(any(mask(:,slice,:) > 0))
             cst{s,7}{slice,2} = contourc(squeeze(mask(:,slice,:)),.5*[1 1]);
        end
    end
    for slice = 1:ct.cubeDim(3)
        if any(any(mask(:,:,slice) > 0))
             cst{s,7}{slice,3} = contourc(squeeze(mask(:,:,slice)),.5*[1 1]);
        end
    end
end
