% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decription: Generate time array formatted for log file creation
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 10/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timeArray = matRad_getTime4log
uhrzeit = datestr(now,'HH:MM:SS');
datum  =  datestr(now,'yyyy-mm-dd');

timeArray = strcat(datum(1:4), datum(6:7), datum(9:10), 'd', uhrzeit(1:2), ...
    uhrzeit(4:5), uhrzeit(7:8), 'h');