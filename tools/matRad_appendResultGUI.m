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

if ~exist('boolOverwrite','var')
    boolOverwrite = false;
end

if ~exist('Identifier','var')
   Identifier = '2';
end

caFieldnames = fieldnames(resultGUItoAppend);

for i = 1:numel(caFieldnames)

   if boolOverwrite
       resultGUI.(caFieldnames{i,1}) =   resultGUItoAppend.(caFieldnames{i,1}); 
   else
       resultGUI.([caFieldnames{i,1} '_' Identifier]) =   resultGUItoAppend.(caFieldnames{i,1}); 
   end
   
end

% group similar fields together
resultGUI = orderfields(resultGUI);

end

