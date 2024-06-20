function [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
% matRad inverse planning wrapper function
%
% call
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln)
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln,wInit)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   wInit:      (optional) custom weights to initialize problems
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
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
if exist('wInit','var')
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln,wInit);
else
    [dij,cst,pln,wInit,optiProb,FLAG_ROB_OPT] = matRad_initOptimization(dij,cst,pln);
end

if ~isfield(pln.propOpt,'optimizer')
    pln.propOpt.optimizer = 'IPOPT';
end

switch pln.propOpt.optimizer
    case 'IPOPT'
        optimizer = matRad_OptimizerIPOPT;
    case 'fmincon'
        optimizer = matRad_OptimizerFmincon;
    otherwise
        warning(['Optimizer ''' pln.propOpt.optimizer ''' not known! Fallback to IPOPT!']);
        optimizer = matRad_OptimizerIPOPT;
end
        
if ~optimizer.IsAvailable()
    matRad_cfg.dispError(['Optimizer ''' pln.propOpt.optimizer ''' not available!']);
end

optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

wOpt = optimizer.wResult;
info = optimizer.resultInfo;

resultGUI = matRad_calcCubes(wOpt,dij);
resultGUI.wUnsequenced = wOpt;
resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;
resultGUI.optiProb = optiProb;

%Robust quantities
if FLAG_ROB_OPT || numel(optiProb.BP.scenarios) > 1
    Cnt = 1;
    for i = find(~cellfun(@isempty,dij.physicalDose))'
        tmpResultGUI = matRad_calcCubes(wOpt,dij,i);
        resultGUI.([pln.bioParam.quantityVis '_' num2str(Cnt,'%d')]) = tmpResultGUI.(pln.bioParam.quantityVis);
        Cnt = Cnt + 1;
    end
end

% unblock mex files
clear mex

end    