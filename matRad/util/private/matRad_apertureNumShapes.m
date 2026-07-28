function n = matRad_apertureNumShapes(apertureInfo, i)
% matRad_apertureNumShapes Number of aperture shapes of beam i.
%
%   Deliberately counts the shape structs rather than reading
%   beam(i).numOfShapes: in VMAT, interpolated (non-DAO) beams carry a
%   computed shape while their numOfShapes (the number of *optimized*
%   shapes) is 0 - and those interpolated apertures should be shown too.
%
% call:
%   n = matRad_apertureNumShapes(apertureInfo, i)
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

n = numel(apertureInfo.beam(i).shape);

end
