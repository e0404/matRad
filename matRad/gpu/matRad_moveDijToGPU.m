function dij = matRad_moveDijToGPU(dij, precision)
% matRad_moveDijToGPU transfers dij dose influence matrix data to the GPU.
%   Moves the dose influence quantities (physicalDose, mAlphaDose,
%   mSqrtBetaDose, mLETDose) to GPU memory. Optionally casts the data to
%   the requested numeric precision before uploading.
%
% call:
%   dij = matRad_moveDijToGPU(dij)
%   dij = matRad_moveDijToGPU(dij, precision)
%
% input:
%   dij         matRad dij struct with host arrays in dose influence fields
%   precision   (optional) target numeric precision as string, e.g.
%               'single' or 'double'. If empty or omitted, no cast is
%               performed.
%
% output:
%   dij         matRad dij struct with GPU arrays in dose influence fields
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

quantities = {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose'};

for i = 1:numel(quantities)
    if isfield(dij, quantities{i})
        if ~isempty(precision)
            dij.(quantities{i}) = cellfun(@(x) gpuArray(cast(x, precision)), dij.(quantities{i}), 'UniformOutput', false);
        else
            dij.(quantities{i}) = cellfun(@(x) gpuArray(x), dij.(quantities{i}), 'UniformOutput', false);
        end
    end
end

end
