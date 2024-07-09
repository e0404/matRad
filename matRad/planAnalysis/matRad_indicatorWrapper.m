function [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI,refGy,refVol)
% matRad indictor wrapper
% 
% call
%   [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI)
%   [dvh,qi] = matRad_indicatorWrapper(cst,pln,resultGUI,refGy,refVol)
%
% input
%   cst:                  matRad cst struct
%   pln:                  matRad pln struct
%   resultGUI:            matRad resultGUI struct
%   refGy: (optional)     array of dose values used for V_XGy calculation
%                         default is [40 50 60]
%   refVol:(optional)     array of volumes (0-100) used for D_X calculation
%                         default is [2 5 95 98]
%                         NOTE: Call either both or none!
%
% output
%   dvh: matRad dvh result struct
%   qi:  matRad quality indicator result struct
%   graphical display of all results
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the matRad configuration instance
matRad_cfg = MatRad_Config.instance();
% Display a deprecation warning for the current function
matRad_cfg.dispDeprecationWarning('The matRad_indicatorWrapper function will be deprecated soon!\nPlan analysis is now handled by matRad_planAnalysis!');

% Initialize an empty cell array for optional arguments to translate the into key-value pairs
args = {};
% Check if 'refVol' variable exists and add it to the arguments if it does
if exist('refVol', 'var') 
    args{end+1,end+2} = {'refVol',refVol};
end

% Check if 'refGy' variable exists and add it to the arguments if it does
if exist('refGy', 'var')
    args{end+1,end+2} = {'refGy',refGy};
end

% Initialize empty structures for ct and stf, required for matRad_planAnalysis
ct = struct();
stf = struct();

% Call the matRad_planAnalysis function with the prepared arguments
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln,args{:});

%Get the return arguments
dvh = resultGUI.dvh;
qi = resultGUI.qi;

end



