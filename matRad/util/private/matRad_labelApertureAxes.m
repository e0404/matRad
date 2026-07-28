function matRad_labelApertureAxes(hAx, leafCoord)
% matRad_labelApertureAxes Axis labels of an aperture plot.
%
% call:
%   matRad_labelApertureAxes(hAx, leafCoord)
%
% input:
%   hAx:       axes to label
%   leafCoord: 'leafNum' or 'physical', see matRad_visApertureInfo
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlabel(hAx, 'horiz. pos. [mm]');

if strcmp(leafCoord, 'physical')
    ylabel(hAx, 'vert. pos. [mm]');
else
    ylabel(hAx, 'leaf pair #');
end

end
