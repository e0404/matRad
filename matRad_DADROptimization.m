function [resultGUI,optimizer] = matRad_DADROptimization(dij,cst,pln,wTilda0)
% matRad inverse planning wrapper function
% 
% call
%   [resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   FLASHindex: array with number of the lines in the cst which should be
%   protected
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

factor = 1;

% issue warning if biological optimization impossible
if sum(strcmp(pln.propOpt.bioOptimization,{'LEMIV_effect','LEMIV_RBExD'}))>0 && (~isfield(dij,'mAlphaDose') || ~isfield(dij,'mSqrtBetaDose')) && strcmp(pln.radiationMode,'carbon')
    warndlg('Alpha and beta matrices for effect based and RBE optimization not available - physical optimization is carried out instead.');
    pln.propOpt.bioOptimization = 'none';
end

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);


% check & adjust objectives and constraints internally for fractionation 
for i = 1:size(cst,1)
    %Compatibility Layer for old objective format
    if isstruct(cst{i,6})
        % DO I HAVE TO ADJUST THIS?
        cst{i,6} = arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6},'UniformOutput',false);
    end
    for j = 1:numel(cst{i,6})
        %In case it is a default saved struct, convert to object
        %Also intrinsically checks that we have a valid optimization
        %objective or constraint function in the end
        obj = cst{i,6}{j};

        % NOT SURE WHAT THIS THUS. STILL NEED TO CHECK
        if ~isa(obj,'matRad_DoseOptimizationFunction') && ~isa(obj,'matRad_DADROptimizationFunction')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
            end
        end
        
        % Checking dose constraints
        if isa(obj,'matRad_DoseOptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
            cst{i,6}{j} = obj;
        elseif  isa(obj,'matRad_DADROptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters());
            cst{i,6}{j} = obj;
            
        end      
        
        
    end
    
    
    
             
   
end

% resizing cst to dose cube resolution 
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                                 dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
dadrTarget = [];

ixTarget   = [];

for i = 1:size(cst,1)
    if ((isequal(cst{i,3},'TARGET') || isequal(cst{i,3},'OAR')) && ~isempty(cst{i,6}))
        V = [V;cst{i,4}{1}];
        
        %Iterate through objectives/constraints
        fDoses = [];
        fDADR = [];
        for fObjCell = cst{i,6}
            if isa(fObjCell{1},'matRad_DoseOptimizationFunction')
                dParams = fObjCell{1}.getDoseParameters();
                %Don't care for Inf constraints
                dParams = dParams(isfinite(dParams));
                %Add do dose list
                
                fDoses = [fDoses dParams];
            elseif isa(fObjCell{1},'matRad_DADROptimizationFunction')
                dadrParams = fObjCell{1}.getDoseParameters();
                
                dadrParams = dadrParams(isfinite(dadrParams));
                
                
                fDADR = [fDADR dadrParams];
                
            end
        end
    
        
        doseTarget = [doseTarget fDoses];
        dadrTarget = [dadrTarget fDADR];

        ixTarget   = [ixTarget i*ones(1,length(fDoses))];
    end
end

%% I don't think this lines are important ( specially the one for dadr)
[doseTarget,i] = max(doseTarget);
[dadrTarget,j] = max(dadrTarget);


ixTarget       = ixTarget(i);
wOnes          = ones(dij.totalNumOfBixels,1);

% modified settings for photon dao
if pln.propOpt.runDAO && strcmp(pln.radiationMode,'photons')
   options.ipopt.max_iter = 50;
   options.ipopt.acceptable_obj_change_tol     = 7e-3; % (Acc6), Solved To Acceptable Level if (Acc1),...,(Acc6) fullfiled

end
% calculate initial beam intensities wInit
if nargin < 4 || isempty(wTilda0)

    if  strcmp(pln.propOpt.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')

        % check if a constant RBE is defined - if not use 1.1
        if ~isfield(dij,'RBE')
            dij.RBE = 1.1;
        end
        bixelWeight =  (doseTarget)/(dij.RBE * mean(dij.physicalDose{1}(V,:)*wOnes)); 
        wInit       = wOnes * bixelWeight; 
        IInit       = [wOnes* 1];

    elseif (strcmp(pln.propOpt.bioOptimization,'LEMIV_effect') || strcmp(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                                    && strcmp(pln.radiationMode,'carbon')

        % retrieve photon LQM parameter
        [ax,bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1);

        if ~isequal(dij.ax(dij.ax~=0),ax(dij.ax~=0)) || ...
           ~isequal(dij.bx(dij.bx~=0),bx(dij.bx~=0))
             matRad_cfg.dispError('Inconsistent biological parameter - please recalculate dose influence matrix!\n');
        end

        for i = 1:size(cst,1)

            for j = 1:size(cst{i,6},2)
                % check if prescribed doses are in a valid domain
                if any(cst{i,6}{j}.getDoseParameters() > 5) && isequal(cst{i,3},'TARGET')
                    matRad_cfg.dispError('Reference dose > 10 Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.\n');
                end

            end
        end

        dij.ixDose  = dij.bx~=0; 

        if isequal(pln.propOpt.bioOptimization,'LEMIV_effect')

               effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
               p            = (sum(dij.mAlphaDose{1}(V,:)*wOnes)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
               q            = -(effectTarget * length(V)) / (sum((dij.mSqrtBetaDose{1}(V,:) * wOnes).^2));
               wInit        = -(p/2) + sqrt((p^2)/4 -q) * wOnes;
               IInit       = [wOnes* 1];

        elseif isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')

               %pre-calculations
               dij.gamma              = zeros(dij.doseGrid.numOfVoxels,1);   
               dij.gamma(dij.ixDose) = dij.ax(dij.ixDose)./(2*dij.bx(dij.ixDose)); 

               % calculate current in target
               CurrEffectTarget = (dij.mAlphaDose{1}(V,:)*wOnes + (dij.mSqrtBetaDose{1}(V,:)*wOnes).^2);
               % ensure a underestimated biological effective dose 
               TolEstBio        = 1.2;
               % calculate maximal RBE in target
               maxCurrRBE = max(-cst{ixTarget,5}.alphaX + sqrt(cst{ixTarget,5}.alphaX^2 + ...
                            4*cst{ixTarget,5}.betaX.*CurrEffectTarget)./(2*cst{ixTarget,5}.betaX*(dij.physicalDose{1}(V,:)*wOnes)));
               wInit    =  ((doseTarget)/(TolEstBio*maxCurrRBE*max(dij.physicalDose{1}(V,:)*wOnes)))* wOnes;
               IInit       = [wOnes* 1];
        end

    else
    %      % retrieve photon LQM parameter
    %     [ax,bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1);
         % I DON'T THINK THIS IS MAKING SENSE


        if ~isempty(doseTarget)
                bixelWeight =  (doseTarget)/(mean(dij.physicalDose{1}(V,:)*wOnes));
        else
            bixelWeight = 15/(mean(dij.physicalDose{1}(V,:)*wOnes)); % this doesn't make much sense, just to try
        end

    %     bixelWeight =  (doseTarget)/(mean(dij.physicalDose{1}(V,:)*wOnes));
        wInit       = [wOnes * bixelWeight];

        dij.ixDose  = dij.physicalDose{1}*wInit~=0; 

        IInit       = [wOnes* 1]; %For this case in the collapsed script, min 7.3293e+11 particles
        pln.propOpt.bioOptimization = 'none';
    end
    
%     IInit(3:end) = IInit(3:end)*10^-2;
    IInit = [IInit*1e3]; %We use 1e3 here because it is 10^9 (standard particles/s) / 10^6 (noramlization of dij / w in particles)
    wTildaInit       = [wInit; IInit];
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% When the precomputed wInit is included as an argument %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else
    
    
    if  strcmp(pln.propOpt.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')

        % check if a constant RBE is defined - if not use 1.1
        if ~isfield(dij,'RBE')
            dij.RBE = 1.1;
        end
        wTildaInit       = wTilda0;

        elseif (strcmp(pln.propOpt.bioOptimization,'LEMIV_effect') || strcmp(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                                    && strcmp(pln.radiationMode,'carbon')

            % retrieve photon LQM parameter
            [ax,bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1);

            if ~isequal(dij.ax(dij.ax~=0),ax(dij.ax~=0)) || ...
               ~isequal(dij.bx(dij.bx~=0),bx(dij.bx~=0))
                 matRad_cfg.dispError('Inconsistent biological parameter - please recalculate dose influence matrix!\n');
            end

            for i = 1:size(cst,1)

                for j = 1:size(cst{i,6},2)
                    % check if prescribed doses are in a valid domain
                    if any(cst{i,6}{j}.getDoseParameters() > 5) && isequal(cst{i,3},'TARGET')
                        matRad_cfg.dispError('Reference dose > 10 Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.\n');
                    end

                end
            end

                dij.ixDose  = dij.bx~=0; 

                if isequal(pln.propOpt.bioOptimization,'LEMIV_effect')
                    
                   wTildaInit       = wTilda0;

                elseif isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')
                   %pre-calculations
                   dij.gamma              = zeros(dij.doseGrid.numOfVoxels,1);   
                   dij.gamma(dij.ixDose) = dij.ax(dij.ixDose)./(2*dij.bx(dij.ixDose)); 
                   wTildaInit       = wTilda0;
                end

        else 

            wTildaInit       = wTilda0;
            wFinal = round(length(wTildaInit)/2);
            wInit = wTildaInit(1:wFinal);
            dij.ixDose  = dij.physicalDose{1}*wInit~=0;
            pln.propOpt.bioOptimization = 'none';
            
    end
end

%Select Projection

switch pln.propOpt.bioOptimization
    case 'LEMIV_effect'
        backProjection = matRad_EffectProjection;
    case 'const_RBExD'
        backProjection = matRad_ConstantRBEProjection;
    case 'LEMIV_RBExD'
        backProjection = matRad_VariableRBEProjection;
    % The default optimization tool will be the FLASH tool    
    case 'none'
        backProjection = matRad_DoseProjection;
    otherwise
        warning(['Did not recognize bioloigcal setting ''' pln.probOpt.bioOptimization '''!\nUsing physical dose optimization!']);
        backProjection = matRad_DoseProjection;
end
        
optiProb = matRad_OptimizationProblemDADR(backProjection,matRad_DADRProjection);

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
        
%optimizer = matRad_OptimizerFmincon;

optimizer = optimizer.optimize(wTildaInit,optiProb,dij,cst);

wFinal = size(dij.physicalDose{1},2);
wTildaOpt = optimizer.wResult;
wOpt = wTildaOpt(1:wFinal);
IOpt = wTildaOpt(1+wFinal:end);
IOpt = IOpt*factor; % Converting from particles/ms to particles/s

% Change the intensity of the IMPT PB's
M = mode(IOpt);
change =(IOpt==M);
IOpt(change) = 1;

% IMPT = find(IOpt==wTildaInit(wFinal+1:end)); %to find the pencil beams whose dose-rate is not changed with the optimizer
%                            % since there's no dose-rate constraint on them
% 
% IOpt(IMPT) = 10^(-2); % I am saying all PBs have 10^8 particle/s intensity

% % Only keep the w's that are greater then 10^6 particles
% [a,aa]=find(wOpt>1);
% 
% wOpt = wOpt(a);
% IOpt = IOpt(a);


info = optimizer.resultInfo;

%% Update of dij to only keep the wanted w's
% dij.beamNum = dij.beamNum(a);
% dij.physicalDose{1} = dij.physicalDose{1};


% dij.physicalDose{1}(:,1:2) = dij.physicalDose{1}(:,1:2)./[10^4    10^4;];
% wOpt(1:2) = wOpt(1:2).*([2.9967e+05    3.1132e+05])';
% IOpt(1:2) = IOpt(1:2).*([ 2.9967e+05    3.1132e+05])';

% resultGUI = matRad_calcCubes(wOpt,dij);
resultGUI = matRad_calcCubes(wOpt,dij);
resultGUI = matRad_calcDADR(dij,resultGUI,IOpt);
resultGUI.I = IOpt;

resultGUI.wUnsequenced = wOpt;
resultGUI.IUnsequenced = IOpt;

resultGUI.usedOptimizer = optimizer;
resultGUI.info = info;

% unblock mex files
clear mex
