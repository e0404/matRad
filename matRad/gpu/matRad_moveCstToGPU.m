function cst = matRad_moveCstToGPU(cst, indexType)
% matRad_moveCstToGPU transfers cst dose influence data to the GPU.
%   Moves all cell arrays stored in column 4 of the cst (dose influence
%   data) to GPU memory. Optionally casts the data to the requested numeric
%   precision before uploading.
%
% call
%   cst = matRad_moveCstToGPU(cst)
%   cst = matRad_moveCstToGPU(cst, precision)
%
% input
%   cst         matRad cst cell array with host arrays in column 4
%   indexType   (optional) target index type as integer, e.g.
%               'int32' or 'uint64'. If empty or omitted, no cast is
%               performed.
%
% output
%   cst         matRad cst cell array with GPU arrays in column 4
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
    indexType = [];
end

for i = 1:size(cst, 1)

    if ~isempty(indexType)
        cst{i, 4} = cellfun(@(x) cast(x, indexType), cst{i, 4}, 'UniformOutput', false);
    end
    cst{i, 4} = cellfun(@gpuArray, cst{i, 4}, 'UniformOutput', false);

end

end
