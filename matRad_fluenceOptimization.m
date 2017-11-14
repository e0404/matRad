function [resultGUI,info] = matRad_fluenceOptimization(dij,cst,pln,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad inverse planning wrapper function
% 
% call
%   [resultGUI,info] = matRad_fluenceOptimization(dij,cst,pln)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   param:      (optional) structure defining additional parameter            
%               e.g. param.logLevel defines the log level
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   info:       struct containing information about optimization
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if exist('param','var')
   if ~isfield(param,'logLevel')
      param.logLevel = 1;
   end
else
   param.logLevel = 1;
end

if ~isdeployed % only if _not_ running as standalone
    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'optimization'))
    
    if param.logLevel == 1
       % get handle to Matlab command window
       mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
       cw          = mde.getClient('Command Window');
       xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
       h_cw        = handle(xCmdWndView,'CallbackProperties');

       % set Key Pressed Callback of Matlab command window
       set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);
    end
end

% initialize global variables for optimizer
global matRad_global_x;
global matRad_global_d;
global matRad_STRG_C_Pressed;
global matRad_objective_function_value;

matRad_global_x                 = NaN * ones(dij.totalNumOfBixels,1);
matRad_global_d                 = NaN * ones(dij.numOfVoxels,1);
matRad_STRG_C_Pressed           = false;
matRad_objective_function_value = [];
  
% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

% adjust objectives and constraints internally for fractionation 
for i = 1:size(cst,1)
    for j = 1:size(cst{i,6},1)
       cst{i,6}(j).dose = cst{i,6}(j).dose/pln.numOfFractions;
    end
end

% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
ixTarget   = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}{1}];
        doseTarget = [doseTarget cst{i,6}.dose];
        ixTarget   = [ixTarget i*ones(1,length([cst{i,6}.dose]))];
    end
end
[doseTarget,i] = max(doseTarget);
ixTarget       = ixTarget(i);
wOnes          = ones(dij.totalNumOfBixels,1);

% set the IPOPT options.
matRad_ipoptOptions;

% modified settings for photon dao
if pln.runDAO && strcmp(pln.radiationMode,'photons')
%    options.ipopt.max_iter = 50;
%    options.ipopt.acceptable_obj_change_tol     = 7e-3; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled
end

% set bounds on optimization variables
options.lb       = zeros(1,dij.totalNumOfBixels);        % Lower bound on the variables.
options.ub       = inf * ones(1,dij.totalNumOfBixels);   % Upper bound on the variables.
funcs.iterfunc   = @(iter,objective,paramter) matRad_IpoptIterFunc(iter,objective,paramter,options.ipopt.max_iter,param);
    
% calculate initial beam intensities wInit
if  strcmp(pln.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
    % check if a constant RBE is defined - if not use 1.1
    if ~isfield(dij,'RBE')
        dij.RBE = 1.1;
    end
    bixelWeight =  (doseTarget)/(dij.RBE * mean(dij.physicalDose{1}(V,:)*wOnes)); 
    wInit       = wOnes * bixelWeight;
        
elseif pln.bioParam.bioOpt
    
    for i = 1:size(cst,1)

        for j = 1:size(cst{i,6},2)
            % check if prescribed doses are in a valid domain
            if cst{i,6}(j).dose > 5 && isequal(cst{i,3},'TARGET')
                matRad_dispToConsole('Reference dose > 5Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.',param,'error');
            end
            
        end
    end
    
    dij.ixDose  = dij.betaX~=0; 
        
    if isequal(pln.bioParam.quantityOpt,'effect')
        
           effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
           p            = (sum(dij.mAlphaDose{1}(V,:)*wOnes)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
           q            = -(effectTarget * length(V)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
           wInit        = -(p/2) + sqrt((p^2)/4 -q) * wOnes;

    elseif isequal(pln.bioParam.quantityOpt,'RBExD')
        
           %pre-calculations
           dij.gamma              = zeros(dij.numOfVoxels,1);   
           dij.gamma(dij.ixDose) = dij.alphaX(dij.ixDose)./(2*dij.betaX(dij.ixDose)); 
            
           % calculate current in target
           CurrEffectTarget = (dij.mAlphaDose{1}(V,:)*wOnes + (dij.mSqrtBetaDose{1}(V,:)*wOnes).^2);
           % ensure a underestimated biological effective dose 
           TolEstBio        = 1.2;
           % calculate maximal RBE in target
           maxCurrRBE = max(-cst{ixTarget,5}.alphaX + sqrt(cst{ixTarget,5}.alphaX^2 + ...
                        4*cst{ixTarget,5}.betaX.*CurrEffectTarget)./(2*cst{ixTarget,5}.betaX*(dij.physicalDose{1}(V,:)*wOnes)));
           wInit    =  ((doseTarget)/(TolEstBio*maxCurrRBE*max(dij.physicalDose{1}(V,:)*wOnes)))* wOnes;
    end
    
else 
    bixelWeight =  (doseTarget)/(mean(dij.physicalDose{1}(V,:)*wOnes)); 
    wInit       = wOnes * bixelWeight;
end

% set optimization options
options.ixForOpt     = find(~cellfun(@isempty, dij.physicalDose))'; 
options.numOfScen    = pln.multScen.numOfScen;
options.scenProb     = pln.multScen.scenProb;
options.bioOpt       = pln.bioParam.bioOpt;
options.quantityOpt  = pln.bioParam.quantityOpt;
options.model        = pln.bioParam.model;

% set callback functions.
[options.cl,options.cu] = matRad_getConstBoundsWrapper(cst,options);   
funcs.objective         = @(x) matRad_objFuncWrapper(x,dij,cst,options);
funcs.constraints       = @(x) matRad_constFuncWrapper(x,dij,cst,options);
funcs.gradient          = @(x) matRad_gradFuncWrapper(x,dij,cst,options);
funcs.jacobian          = @(x) matRad_jacobFuncWrapper(x,dij,cst,options);
funcs.jacobianstructure = @( ) matRad_getJacobStruct(dij,cst);

% Run IPOPT.
[wOpt, info]            = ipopt(wInit,funcs,options);

% calc dose and reshape from 1D vector to 2D array
matRad_dispToConsole('Calculating final cubes...\n',param,'info');
resultGUI = matRad_calcCubes(wOpt,dij,cst);
resultGUI.wUnsequenced = wOpt;

% calc individual scenarios
if pln.multScen.numOfScen > 1
   Cnt = 1;
   for i = options.ixForOpt
      tmpResultGUI = matRad_calcCubes(wOpt,dij,cst,i);
      resultGUI.([pln.bioParam.quantityVis '_' num2str(Cnt,'%d')]) = tmpResultGUI.(pln.bioParam.quantityVis);
      Cnt = Cnt + 1;
   end      
end

% unset Key Pressed Callback of Matlab command window
if ~isdeployed && param.logLevel == 1
    set(h_cw, 'KeyPressedCallback',' ');
end

% clear global variables
clearvars -global matRad_global_x matRad_global_d matRad_objective_function_value matRad_STRG_C_Pressed;

% unblock mex files
clear mex