function tf = matRad_isgpuarrayCompat(x)
% matRad wrapper around Matlab's isgpuarray function for compatibility.
%   isgpuarray belongs to the Parallel Computing Toolbox, so calling it
%   fails in Matlab installations without the toolbox and in Octave. This
%   wrapper checks the class of the input instead, which works everywhere
%   and is false for any input that cannot be a gpuArray.
%
% call:
%   tf = matRad_isgpuarrayCompat(x)
%
% input:
%   x       any variable
%
% output:
%   tf      true if x is a gpuArray, false otherwise
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

tf = isa(x, 'gpuArray');

end
