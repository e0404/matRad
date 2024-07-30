function result = matRad_ispropCompat(obj,prop)
% matRad function mimicking Matlab's properties for compatibility with
% Octave 6 in classdef files (avoids a parse error in the file)
%
% call
%   result = matRad_ispropCompat(obj)
%
% input
%   obj         object (classdef) to check for property
%   prop        property to check for
%
% output
%   result      true if property exists, false otherwise
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team. 
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

if matRad_cfg.isOctave && str2double(matRad_cfg.envVersion(1)) <= 6
    result = any(strcmp(matRad_getPropsCompat(obj),prop));
else
    result = isprop(obj,prop);
end

end