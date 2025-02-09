function writeRunFile(~, fName)
% FRED helper to write file for simulation
% call
%   writeRunFile(fName)
% 
% input
%   fName: string specifying the file path and name for saving the data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

fID = fopen(fName, 'w');

try
    % Include regions and plan delivery routine
    fprintf(fID, 'include: inp/regions/regions.inp\n');
    fprintf(fID, 'include: inp/plan/planDelivery.inp\n');
catch
    matRad_cfg.dispError('Failed to write run file');
end

fclose(fID);
end