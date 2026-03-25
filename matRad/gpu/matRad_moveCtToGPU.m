function ct = matRad_moveCtToGPU(ct, precision)
% matRad_moveCtToGPU transfers ct cube data to the GPU.
%   Moves ct.cubeHU and ct.cube cell arrays to GPU memory. Optionally
%   casts the data to the requested numeric precision before uploading.
%
% call:
%   ct = matRad_moveCtToGPU(ct)
%   ct = matRad_moveCtToGPU(ct, precision)
%
% input:
%   ct          matRad ct struct with host arrays in cubeHU and/or cube
%   precision   (optional) target numeric precision as string, e.g.
%               'single' or 'double'. If empty or omitted, no cast is
%               performed.
%
% output:
%   ct          matRad ct struct with GPU arrays in cubeHU and/or cube
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

if ~isempty(precision)
    if isfield(ct, 'cubeHU')
        ct.cubeHU = cellfun(@(x) cast(x, precision), ct.cubeHU, 'UniformOutput', false);
    end
    if isfield(ct, 'cube')
        ct.cube = cellfun(@(x) cast(x, precision), ct.cube, 'UniformOutput', false);
    end
end

if isfield(ct, 'cubeHU')
    ct.cubeHU = cellfun(@gpuArray, ct.cubeHU, 'UniformOutput', false);
end

if isfield(ct, 'cube')
    ct.cube = cellfun(@gpuArray, ct.cube, 'UniformOutput', false);
end

end
