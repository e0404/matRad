function [resultGUI, result] = matRad_accumulateCubesMixMod(resultGUI, pln,ct)
% function to merge Mixed modality resultGUI cells into one resultGUI for
% visualisation
% 
% call
%   [resultGUI, result] = matRad_accumulateCubesMixMod(resultGUI, pln,ct)
%
% input
%   resultGUI:              matRads resultGUI Cell array from Mixed
%                           modality optimization
%   pln:                    matRad pln stuct
%   ct:                     matRad ct structure 
%                           
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
matRad_cfg = MatRad_Config.instance();
if ~iscell(resultGUI)
    matRad_cfg.dispWarning('resultGUI is already a structure');
    return;
end
if numel(pln) ~= numel(resultGUI)
    matRad_cfg.dispWarning('dimension mismatch between pln and resultGUI ');
end

quantityOpt = pln(1).bioParam.quantityOpt;
result = resultGUI;
resultGUI = [];
resultGUI.totalQO = zeros(ct.cubeDim);
for mod = 1 : size(result,2)
        modName = ['mod', num2str(mod)];
        for st = 1 : size(result,1)
            stName = ['st',num2str(st)];
            fNames = fieldnames(result{st,mod});
            for i = 1: numel(fNames)
                tn = [modName,stName,fNames{i}];
                resultGUI.(tn) = result{st,mod}.(fNames{i});
            end
        end
         tn = [modName quantityOpt];
        resultGUI.(tn) = result{st,mod}.(quantityOpt) * pln(mod).numOfFractions;
        resultGUI.totalQO = resultGUI.totalQO + resultGUI.(tn);
end


end