function resultGUI = matRad_calcCubes(w,dij,scenNum)
% matRad computation of all cubes for the resultGUI struct
% which is used as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubes(w,dij)
%   resultGUI = matRad_calcCubes(w,dij,scenNum)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   scenNum: optional: number of scenario to calculated (default 1)
%
% output
%   resultGUI: matRad result struct
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
if nargin < 3
    scenNum = 1;
end

% First calculate on dose grid
resultGUI = matRad_calcCubesDoseGrid(w,dij,scenNum);

% interpolation if dose grid does not match ct grid
if isfield(dij,'ctGrid') && any(dij.ctGrid.dimensions~=dij.doseGrid.dimensions)
    myFields = fieldnames(resultGUI);
    for i = 1:numel(myFields)
        if numel(resultGUI.(myFields{i})) == dij.doseGrid.numOfVoxels

            % interpolate!
            resultGUI.(myFields{i}) = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                resultGUI.(myFields{i}), ...
                dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);

        end
    end
end

end

