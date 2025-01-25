function gridStruct = matRad_getWorldAxes(gridStruct)
% matRad function to compute and store world coordinates into ct.x
%
% call
%   gridStruct = matRad_getWorldAxes(gridStruct)
%
%   gridStruct:         can be ct, dij.doseGrid,dij.ctGrid
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

update = false;
if nargin>0
    if ~ (isfield(gridStruct,{'x','y','z'}))
        update = true;
    else
        if isempty(gridStruct.x)||isempty(gridStruct.y)||isempty(gridStruct.z)
            update = true;
        end
    end

    if update
        %
        if isfield(gridStruct,'cubeDim')
            gridStruct.dimensions = gridStruct.cubeDim;
        end
        % check if dicominfo exists
        if isfield(gridStruct, 'dicomInfo') && isfield(gridStruct.dicomInfo,'ImagePositionPatient')
            firstVox = gridStruct.dicomInfo.ImagePositionPatient;
        else
            firstVox = - (gridStruct.dimensions([2 1 3])./2).*[gridStruct.resolution.x gridStruct.resolution.y gridStruct.resolution.z] ;
        end

        gridStruct.x = firstVox(1) + gridStruct.resolution.x*(0:gridStruct.dimensions(2)-1);
        gridStruct.y = firstVox(2) + gridStruct.resolution.y*(0:gridStruct.dimensions(1)-1);
        gridStruct.z = firstVox(3) + gridStruct.resolution.z*(0:gridStruct.dimensions(3)-1);
    end
else
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Cannot compute world coordinates without matRad ct or grid structure');
end

end
