function baseData = matRad_checkBaseData(baseData)
% script to switch from APM (machine.data.Z is a struct) to "original"
% calculation (machine.data.Z is double)
%
% call
%   baseData = matRad_checkBaseData(baseData)
%
% input
%   baseData:       matRad machine base data
%
% output
%   baseData        updated base data with new depths if possible
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(baseData(1).Z) && isfield(baseData(1).Z,'profileORG')
    for i = 1:length(baseData)
        baseData(i).Z = baseData(i).Z.profileORG;
    end
elseif isstruct(baseData(1).Z) && isfield(baseData(1).Z,'doseORG')
    for i = 1:length(baseData)
        baseData(i).Z = baseData(i).Z.doseORG;
    end
elseif isstruct(baseData(1).Z) && (~isfield(baseData(1).Z,'profileORG') && ~isfield(baseData(1).Z,'doseORG'))
    warning('No original depths available in base data. Nothing changed.');
else
    warning('Base data depths are already in the desired format.');
end
end