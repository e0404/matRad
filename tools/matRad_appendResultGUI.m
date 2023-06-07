function resultGUI = matRad_appendResultGUI(resultGUI,resultGUItoAppend,boolOverwrite,Identifier)
% function to merge two seperate resultGUI structs into one for
% visualisation
% 
% call
%   resultGUI = matRad_mergeResultGUIs(resultGUI,resultGUIrob)
%
% input
%   resultGUI:              matRads resultGUI struct
%   resultGUItoAppend:      resultGUI struct which will be appendet
%   boolOverwrite:          if true existing fields be overwritten in case
%                           they already exist
%
% output
%   resultGUI:        matRads resultGUI struct
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Identifier','var')
    Identifier = '';
else
    Identifier = ['_' Identifier];
end

if ~exist('boolOverwrite','var')
    boolOverwrite = false;
end

matRad_cfg = MatRad_Config.instance();

caFieldnames = fieldnames(resultGUItoAppend);

for i = 1:numel(caFieldnames)

    currFieldName = caFieldnames{i,1};
    fullFieldName = [caFieldnames{i,1} Identifier];

    if boolOverwrite || ~isfield(resultGUI,fullFieldName)
        resultGUI.(fullFieldName) =  resultGUItoAppend.(currFieldName);        
    else
        matRad_cfg.dispWarning('Field ''%s'' exists and overwriting was disabled. Results will not be appended to resultGUI!',fullFieldName);
    end
end

% group similar fields together
resultGUI = orderfields(resultGUI);

end

