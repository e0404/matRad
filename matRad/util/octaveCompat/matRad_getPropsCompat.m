function p = matRad_getPropsCompat(obj)
% matRad function mimicking Matlab's properties for compatibility with
% Octave 6 in classdef files (avoids a parse error in the file)
%
% call
%   p =  matRad_getPropsCompat(obj)
%
% input
%   obj         object (classdef) to get properties from
%
% output
%   p           properties of the object
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
try
    p = properties(obj);
catch ME
    p = [];
end
end