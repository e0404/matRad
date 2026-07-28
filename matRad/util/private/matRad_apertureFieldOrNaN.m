function val = matRad_apertureFieldOrNaN(s, fieldName)
% matRad_apertureFieldOrNaN Scalar field of an aperture info struct, or NaN.
%
%   Aperture info structs are grown by the sequencer and the optimizer, so
%   which fields exist depends on how far a plan has been taken - the gantry
%   angle and the delivery timing, for instance, are VMAT-only. Reading them
%   through this function lets a caller display what is there and leave out
%   what is not.
%
% call:
%   val = matRad_apertureFieldOrNaN(s, fieldName)
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

val = NaN;

if isfield(s, fieldName) && ~isempty(s.(fieldName))
    val = s.(fieldName)(1);
end

end
