function matRad_checkEnvDicomRequirements(env)
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


if strcmp(env,'OCTAVE')
    try
        pkg load dicom;
        pkg load image;
    catch
        error('Octave needs dicom and image packages from Octave Forge for dicom import!');
    end
else
    if ~license('checkout','image_toolbox')
        error('image processing toolbox and/or corresponding licence not available');
    end      
end

end

