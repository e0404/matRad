function coord = matRad_world2cubeCoords(wCoord, gridStruct)
% matRad function to convert world coordinates to cube coordinates
% 
% call
%   coord = world2cubeCoords(wCoord, ct)
%
%
%   wCoord:         world coordinates (x,y,z)[mm]
%   gridStruct:     can be matRad ct, dij.doseGrid, or the ctGrid
%                   required fields x,y,x,dimensions,resolution
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
           
        gridStruct = matRad_getWorldAxes(gridStruct);
        
        % world bounds
        boundx = [min(gridStruct.x)-gridStruct.resolution.x/2 max(gridStruct.x)+gridStruct.resolution.x/2];
        boundy = [min(gridStruct.y)-gridStruct.resolution.y/2 max(gridStruct.y)+gridStruct.resolution.y/2];
        boundz = [min(gridStruct.z)-gridStruct.resolution.z/2 max(gridStruct.z)+gridStruct.resolution.z/2];

        % check if queried coordinate is within worldcube
        if (wCoord(1)>=boundx(1) && wCoord(1)<=boundx(2)) && ...
           (wCoord(2)>=boundy(1) && wCoord(2)<=boundy(2)) && ...
           (wCoord(3)>=boundz(1) && wCoord(3)<=boundz(2))
                
            % find the closest world coord index in gridStruct.x/y/z
            % calc absolute differences and locate smallest difference
            [~,coord(1)]=min(abs(gridStruct.x-wCoord(1)));
            [~,coord(2)]=min(abs(gridStruct.y-wCoord(2)));
            [~,coord(3)]=min(abs(gridStruct.z-wCoord(3)));
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
