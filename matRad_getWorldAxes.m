function ct = matRad_getWorldAxes(ct)
% matRad function to compute and store world coordinates int ct.x
% 
% call
%   ct = matRad_getWorldAxes(ct)
%
% References
%   -
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
matRad_cfg = MatRad_Config.instance();

if nargin>0
    %
    % check if dicominfo exists
    if isfield(ct, 'dicomInfo') && isfield(ct.dicomInfo,'ImagePositionPatient')
        firstVox = ct.dicomInfo.ImagePositionPatient;
    else 
        firstVox = - (ct.cubeDim./2).*[ct.resolution.x ct.resolution.y ct.resolution.z] ;
    end 
    
    ct.x = firstVox(1) + ct.resolution.x*[0:ct.cubeDim(2)-1] ;
    ct.y = firstVox(2) + ct.resolution.y*[0:ct.cubeDim(1)-1] ;
    ct.z = firstVox(3) + ct.resolution.z*[0:ct.cubeDim(3)-1] ;
else
    matRad_cfg.dispWarning('Cannot compute world coordinates without matRad ct structure');
end

end