function machine = matRad_loadMachine(pln,filepath)
% matRad_loadMachine load a machine base data file from pln struct
%
% call
%   machine = matRad_loadMachine(pln,filepath)
%   machine = matRad_loadMachine(pln)
%
% input
%   pln:            matRad plan meta information struct
%   filepath:       (optional) filepath where to look for machine. default
%                   is basedata in the matRad root folder 
%
% output
%   machine:        matRad machine struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team. 
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
if isfield(pln, 'radiationMode')
    if isfield(pln, 'machine')
        fileName = [pln.radiationMode '_' pln.machine];
    else
        fileName = [pln.radiationMode '_Generic'];
        matRad_cfg.dispWarning('No machine name given, loading generic machine');
    end
else
    matRad_cfg.dispError('No radiation mode given in pln');
end

if ~exist('filepath','var')
    filepath = [matRad_cfg.matRadRoot filesep 'basedata' filesep fileName];
end

try
    m = load(filepath, 'machine');
    machine = m.machine; % strip first layer of loaded struct for convenience
catch
    matRad_cfg.dispError('Could not find the following machine file: %s\n',fileName);
end
end

