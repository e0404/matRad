function sigmaRashi = matRad_calcSigmaRashi(energy,rangeShifter,SSD)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of additional beam broadening due to the use of range shifters (only for protons)
% 
% call
%   sigmaRashi = matRad_calcSigmaRashi(rangeShifter,SSD)
%
% input
%   energy:       initial particle energy
%   rangeShifter: structure defining range shifter geometry
%   SSD:          source to surface distance
%
% output
%   sigmaRashi:   sigma of range shifter (to be added ^2) in mm
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

warning('Range shifter calculation not included in public matRad release');

sigmaRashi = 0;