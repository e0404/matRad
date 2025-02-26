function machine = matRad_loadMachine(pln)
% matRad_loadMachine load a machine base data file from pln struct. 
%   Looks for the machine file from pln in the basedata folder and in the provided user folders.
%
% call
%   machine = matRad_loadMachine(pln)
%
% input
%   pln:            matRad plan meta information struct 
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

userfolders = cellfun(@(uf) strcat(uf,filesep,'machines',filesep),matRad_cfg.userfolders,'UniformOutput',false);

folders = [{[matRad_cfg.matRadSrcRoot filesep 'basedata' filesep]} userfolders(:)'];

foundData = cellfun(@(folder) exist(fullfile(folder, [fileName '.mat']), 'file'), folders);
foundIx = find(foundData, 1, 'first');

if isempty(foundIx)
    matRad_cfg.dispError('Could not find the following machine file: %s',fileName);
end

filepath = fullfile(folders{foundIx}, [fileName '.mat']);

try
    m = load(filepath, 'machine');
    machine = m.machine; % strip first layer of loaded struct for convenience
catch
    matRad_cfg.dispError('Could not load the following machine file: %s',fileName);
end
end

