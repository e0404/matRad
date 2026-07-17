function available = matRad_checkEnvImageProcessingRequirements(env)
% matRad function to check if image processing functionality is available,
% i.e., the Image Processing Toolbox (MATLAB) or the image package (Octave)
%
% call:
%   available = matRad_checkEnvImageProcessingRequirements()
%   available = matRad_checkEnvImageProcessingRequirements(env)
%
% input:
%   env:         (optional) environment string 'MATLAB' or 'OCTAVE'
%
% output:
%   available:   true if the image processing functionality is available
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

available = true;

matRad_cfg = MatRad_Config.instance();

if nargin < 1
    isOctave = matRad_cfg.isOctave;
else
    isOctave = strcmp(env, 'OCTAVE');
end

if isOctave
    try
        pkg load image;
    catch
        available = false;
    end
else
    if ~license('checkout', 'image_toolbox')
        available = false;
    end
end

end
