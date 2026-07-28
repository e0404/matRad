function c = matRad_apertureWeightColor(weight, wMax)
% matRad_apertureWeightColor Background colour encoding an aperture weight.
%
%   Light grey for an unweighted aperture, saturated red for the heaviest
%   one, so that the relative importance of the apertures of a plan can be
%   read off the plots at a glance.
%
% call:
%   c = matRad_apertureWeightColor(weight, wMax)
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

color = [0.2:0.01:0.8; 0.2:0.01:0.8; 0.2:0.01:0.8]';
color = flipud(color);
color(:, 3) = 0;
color(:, 2) = 0;

colorInd = max(ceil((weight / wMax) * (size(color, 1) - 1) + eps), 1);
colorInd = min(colorInd, size(color, 1));
c = color(colorInd, :);

end
