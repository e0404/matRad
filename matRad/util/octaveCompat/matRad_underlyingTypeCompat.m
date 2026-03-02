function utype = matRad_underlyingTypeCompat(x)
% matRad function to obtain the type of a numeric datatype.
%   Matlab has the function underlyingType to robustly determine the type
%   of numerical data, since class might return, for example, "gpuArray" in
%   case of data stored on the GPU compared to "single", or "double" for
%   standard arrays. This wraps the function for Octave compatibility.
%
% call
%   utype = matRad_underlyingTypeCompat(x)
%
% input
%   x       object to check for datatype
%
% output
%   utype   datatype
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
    utype = underlyingType(x);
else
    utype = class(x);
end

end
