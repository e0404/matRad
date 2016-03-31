function [resultGUI,info] = matRad_fluenceOptimization(dij,cst,pln)
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
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isdeployed % only if _not_ running as standalone
    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'optimization'))
    
    % get handle to Matlab command window
    mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
    cw          = mde.getClient('Command Window');
    xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
    h_cw        = handle(xCmdWndView,'CallbackProperties');

    % set Key Pressed Callback of Matlab command window
    set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);

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
        V = [V;cst{i,4}];
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
    options.ipopt.max_iter = 30;
end

% set bounds on optimization variables
options.lb              = zeros(1,dij.totalNumOfBixels);        % Lower bound on the variables.
options.ub              = inf * ones(1,dij.totalNumOfBixels);   % Upper bound on the variables.
funcs.iterfunc          = @(iter,objective,paramter) matRad_IpoptIterFunc(iter,objective,paramter,options.ipopt.max_iter);
    
% calculate initial beam intensities wInit
if strcmp(pln.bioOptimization,'none')
    
    bixelWeight =  (doseTarget)/(mean(dij.physicalDose{1}(V,:)*wOnes)); 
    wInit       = wOnes * bixelWeight;

else (isequal(pln.bioOptimization,'effect') || isequal(pln.bioOptimization,'RBExD')) ... 
        && isequal(pln.radiationMode,'carbon');

    % check if you are running a supported rad
    dij.ax   = zeros(dij.numOfVoxels,1);
    dij.bx   = zeros(dij.numOfVoxels,1);
    
    for i = 1:size(cst,1)
        
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
             dij.ax(cst{i,4}) = cst{i,5}.alphaX;
             dij.bx(cst{i,4}) = cst{i,5}.betaX;
        end
        
        for j = 1:size(cst{i,6},2)
            % check if prescribed doses are in a valid domain
            if cst{i,6}(j).dose > 30 && isequal(cst{i,3},'TARGET')
                error('Reference dose > 30Gy for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.');
            end
            
        end
    end
     
    if isequal(pln.bioOptimization,'effect')
        
           effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
           p            = (sum(dij.mAlphaDose(V,:)*wOnes)) / (sum((dij.mSqrtBetaDose(V,:) * wOnes).^2));
           q            = -(effectTarget * length(V)) / (sum((dij.mSqrtBetaDose(V,:) * wOnes).^2));
           wInit        = -(p/2) + sqrt((p^2)/4 -q) * wOnes;

    elseif isequal(pln.bioOptimization,'RBExD')
        
           %pre-calculations
           dij.gamma      = zeros(dij.numOfVoxels,1);
           idx            = dij.bx~=0; 
           dij.gamma(idx) = dij.ax(idx)./(2*dij.bx(idx)); 
            
           roughRBE = 3;
           wInit    =  (doseTarget)/(roughRBE*mean(dij.physicalDose(V,:)*wOnes)) * wOnes; 
    end
end

% set callback functions.
[options.cl,options.cu] = matRad_getConstBounds(cst,pln.bioOptimization);   
funcs.objective         = @(x) matRad_objFuncWrapper(x,dij,cst,pln.bioOptimization);
funcs.constraints       = @(x) matRad_constFunc(x,dij,cst,pln.bioOptimization);
funcs.gradient          = @(x) matRad_gradFuncWrapper(x,dij,cst,pln.bioOptimization);
funcs.jacobian          = @(x) matRad_jacobFunc(x,dij,cst,pln.bioOptimization);
funcs.jacobianstructure = @( ) matRad_getJacobStruct(dij,cst);

% Run IPOPT.
[wOpt, info]           = ipopt(wInit,funcs,options);
resultGUI.w            = wOpt;
resultGUI.wUnsequenced = wOpt;

% calc dose and reshape from 1D vector to 2D array
resultGUI.physicalDose = reshape(dij.physicalDose{1}*resultGUI.w,dij.dimensions);

if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
                            && strcmp(pln.radiationMode,'carbon')

    fprintf('Calculating alpha/beta/effect/cube...');

    a_x = zeros(size(resultGUI.physicalDose));
    b_x = zeros(size(resultGUI.physicalDose));
    for  i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,4}) = cst{i,5}.alphaX;
            b_x(cst{i,4}) = cst{i,5}.betaX;
        end
    end
    
    resultGUI.effect = (dij.mAlphaDose*resultGUI.w+(dij.mSqrtBetaDose*resultGUI.w).^2);
    resultGUI.effect = reshape(resultGUI.effect,dij.dimensions);
    
    resultGUI.RBExDose     = zeros(size(resultGUI.effect));
    ix                     = resultGUI.effect>0;
    resultGUI.RBExDose(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.effect(ix)) - a_x(ix))./(2.*b_x(ix)));
    resultGUI.RBE          = resultGUI.RBExDose./resultGUI.physicalDose;
   
    AlphaDoseCube    = dij.mAlphaDose * resultGUI.w;
    resultGUI.alpha  = (reshape(AlphaDoseCube,dij.dimensions))./resultGUI.physicalDose;
    SqrtBetaDoseCube = dij.mSqrtBetaDose * resultGUI.w;
    resultGUI.beta   = ((reshape(SqrtBetaDoseCube,dij.dimensions))./resultGUI.physicalDose).^2;
    
    fprintf(' done!\n');
 
end

% unset Key Pressed Callback of Matlab command window
if ~isdeployed
    set(h_cw, 'KeyPressedCallback',' ');
end

clearvars -global matRad_global_x matRad_global_d matRad_ALT_C_Pressed;

