function x = matRad_gatherCompat(x)
% matRad wrapper around Matlab's gather function for Octave compatibility.
%   In Matlab, gather transfers a gpuArray from device to host memory.
%   Octave has no GPU array support, so this function returns the input
%   unchanged.
%
% call:
%   x = matRad_gatherCompat(x)
%
% input:
%   x       array (gpuArray on Matlab, regular array on Octave)
%
% output:
%   x       array in host memory
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

matRad_cfg = MatRad_Config.instance();

if matRad_cfg.isMatlab
    x = gather(x);
end

end
