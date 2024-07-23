function coord = matRad_world2cubeCoords(wCoord, ct)
% matRad function to convert world coordinates to cube coordinates
% 
% call
%   coord = world2cubeCoords(wCoord, ct)
%
%
%   wCoord:         world coordinates (x,y,z)[mm]
%   ct    :         matRad ct structure  
%
%   coord :         cube index (x,y,z)     
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

if nargin == 2 && ~isempty(wCoord)
    
        if ~ (isfield(ct,'x') && isfield(ct,'y') && isfield(ct,'z') )
            ct = matRad_getWorldAxes(ct);
        end
        % world bounds
        boundx = [min(ct.x)-ct.resolution.x/2 max(ct.x)+ct.resolution.x/2];
        boundy = [min(ct.y)-ct.resolution.y/2 max(ct.y)+ct.resolution.y/2];
        boundz = [min(ct.z)-ct.resolution.z/2 max(ct.z)+ct.resolution.z/2];

        if (wCoord(1)>=boundx(1) && wCoord(1)<=boundx(2)) && ...
           (wCoord(2)>=boundy(1) && wCoord(2)<=boundy(2)) && ...
           (wCoord(3)>=boundz(1) && wCoord(3)<=boundz(2))
                
            [~,coord(1)]=min(abs(ct.x-wCoord(1)));
            [~,coord(2)]=min(abs(ct.y-wCoord(2)));
            [~,coord(3)]=min(abs(ct.z-wCoord(3)));
        else 
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('Queried world coordinate is outside the ct');
            
        end

%     catch
% 
%     end
else
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Cannot compute cube coordinates without matRad input coordinates or ct structure');
end

end
