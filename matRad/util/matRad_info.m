function message = matRad_info()
% matRad function to get information message
% 
% call
%   message = matRad_info()
%
% input
%
% output
%   message:  An Information message about matRad
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Copyright
currYear = datestr(now,'yyyy');
message = ['Copyright (C) ' currYear ' The matRad developers @ DKFZ Group Radiotherapy Optimization (e0404)'];

%Warranty clause and license
message = [message newline newline ...
    'matRad is free software and NOT A MEDICAL PRODUCT!' ...
    newline 'matRad comes WITHOUT ANY WARRANTY and does not guarantee FITNESS FOR A PARTICULAR PURPOSE.' ...
    newline 'Use matRad ONLY for research and education at your own risk.' ...
    newline 'For license conditions see LICENSE.md. Third-party software used in matRad is subject to their respective licenses.'];             

%Websites
message = [message newline newline ...
    'Check www.matRad.org and github.com/e0404/matRad for more information.'];

end