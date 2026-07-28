function panels = matRad_aperturePanels(apertureInfo, numOfBeams)
% matRad_aperturePanels Flattens all (beam, shape) pairs into one list.
%
%   The list is in delivery order: the beams are kept in the order they are
%   stored in, which for VMAT is the order in which the arc is traversed
%   (sorting by gantry angle would break reverse arcs and multi-arc plans).
%
% call:
%   panels = matRad_aperturePanels(apertureInfo, numOfBeams)
%
% output:
%   panels: numPanels-by-2 array of [beam index, shape index]
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

panels = zeros(0, 2);

for i = 1:numOfBeams
    for j = 1:matRad_apertureNumShapes(apertureInfo, i)
        panels(end + 1, :) = [i j]; %#ok<AGROW>
    end
end

end
