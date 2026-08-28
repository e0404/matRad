function cstBodyIndex = matRad_findBodyStructureCST(cst, bodyStructureName)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB function to find the index of the radiotherapy structure body.
%
% Note: Body structure is the only structure that has to be contoured.
%
% call
%   cstBodyIndex = matRad_findBodyStructureCST(cst)
%
% input
%   cst
%   bodyStructureName
%
% output:
%   cstBodyIndex
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 05/2019
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

matRad_cfg = MatRad_Config.instance();

cstBodyIndex = 1;

while ~strcmpi(cst{cstBodyIndex,2}, bodyStructureName)
    cstBodyIndex = cstBodyIndex +1;
    if cstBodyIndex > size(cst,1)
        matRad_cfg.dispError('No body structure contoured or structure not found! Note: Body structure has to be named BODY (case insensitive).');
    end
end