function coord = world2cubeCoords(wCoord, ct)
% matRad function to convert world coordinates to cube coordinates
% 
% call
%   coord = world2cubeCoords(wCoord, ct)
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
coord = [];
if nargin == 2
    try
        if ~ (isfield(ct,'x') && isfield(ct,'y') && isfield(ct,'z') )
            ct = matRad_getWorldAxes(ct);
        end
        
        [~,coord(2)]=min(abs(ct.y-wCoord(2)));
        [~,coord(1)]=min(abs(ct.x-wCoord(1)));
        [~,coord(3)]=min(abs(ct.z-wCoord(3)));

    catch

    end
end

end


