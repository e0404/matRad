function available = matRad_checkEnvDicomRequirements(env)
% matRad function to check if requirements for dicom import / export are
% given. Throws an error if requirements not met
% 
% call
%   matRad_checkEnvDicomRequirements(env)
%
% input
%   env:         folder to be scanned
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
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
    isOctave = strcmp(env,'OCTAVE');
end

if isOctave
    try
        pkg load dicom;
        pkg load image;
    catch
        available = false;
    end
else
    if ~license('checkout','image_toolbox')
        available = false;
    end      
end

end

