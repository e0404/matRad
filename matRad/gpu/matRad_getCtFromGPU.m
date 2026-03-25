function ct = matRad_getCtFromGPU(ct, precision)
% matRad_getCtFromGPU transfers ct cube data from GPU to host.
%   Gathers ct.cubeHU and ct.cube cell arrays from GPU memory back to host
%   memory. Optionally casts the data to the requested numeric precision
%   after gathering.
%
% call:
%   ct = matRad_getCtFromGPU(ct)
%   ct = matRad_getCtFromGPU(ct, precision)
%
% input:
%   ct          matRad ct struct with GPU arrays in cubeHU and/or cube
%   precision   (optional) target numeric precision as string, e.g.
%               'single' or 'double'. If empty or omitted, no cast is
%               performed.
%
% output:
%   ct          matRad ct struct with host arrays in cubeHU and/or cube
%
% References
%   -
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

if nargin < 2
    precision = [];
end

if isfield(ct, 'cubeHU')
    ct.cubeHU = cellfun(@matRad_gatherCompat, ct.cubeHU, 'UniformOutput', false);
    if ~isempty(precision)
        ct.cubeHU = cellfun(@(x) cast(x, precision), ct.cubeHU, 'UniformOutput', false);
    end
end

if isfield(ct, 'cube')
    ct.cube = cellfun(@matRad_gatherCompat, ct.cube, 'UniformOutput', false);
    if ~isempty(precision)
        ct.cube = cellfun(@(x) cast(x, precision), ct.cube, 'UniformOutput', false);
    end
end

end
