function machineList = matRad_getAvailableMachines(modalities)
% matRad_loadMachine load a machine base data file from pln struct. 
%   Looks for the machine fileuct.  in the basedata folder and in the provided user folders.
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
% Copyright 2024 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    modalities = {'photons','protons','helium','carbon'};
end

if ischar(modalities)
    modalities = {modalities};
end

matRad_cfg = MatRad_Config.instance();

machineList = containers.Map();

for i = 1:length(modalities)
    pattern = [modalities{i} '_*.mat'];
    machineFiles = dir([matRad_cfg.matRadSrcRoot filesep 'basedata' filesep pattern]);
    for f = 1:length(matRad_cfg.userfolders)
        machineFiles = [machineFiles; dir([matRad_cfg.userfolders{f} filesep 'machines' filesep pattern])];
    end

    if ~isempty(machineFiles)
        machineNames = cell(1,length(machineFiles));
        for j = 1:length(machineFiles)
            machineNames{j} = machineFiles(j).name(numel(modalities{i})+2:end-4);
        end
        machineList(modalities{i}) = machineNames;
    else
        machineList(modalities{i}) = {};
    end
end