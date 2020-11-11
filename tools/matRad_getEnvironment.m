function [env, versionString] = matRad_getEnvironment()
% matRad function to get the software environment matRad is running on
% 
% call
%   [env, versionString] = matRad_getEnvironment()
%
% input
%   -
%   
% output
%   env:            outputs either 'MATLAB' or 'OCTAVE' as string
%   versionString:  returns the version number as string
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reduce runtime if function is called multiple times

matRad_cfg = MatRad_Config.instance();

env = matRad_cfg.env;
versionString = matRad_cfg.envVersion;

end


