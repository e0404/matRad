function [resultGUI,optimizer] = matRad_fluenceOptimizationJO(dij,cst,pln,wInit)
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

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

% check & adjust objectives and constraints internally for fractionation 
for i = 1:size(cst,1)
    %Compatibility Layer for old objective format
    if isstruct(cst{i,6})
        cst{i,6} = arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6},'UniformOutput',false);
    end
    for j = 1:numel(cst{i,6})
        
        obj = cst{i,6}{j};        
        
        %In case it is a default saved struct, convert to object
        %Also intrinsically checks that we have a valid optimization
        %objective or constraint function in the end
        if ~isa(obj,'matRad_DoseOptimizationFunction')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
            end
        end

        if ~strcmp(pln.radiationMode, 'MixMod')
           obj = obj.setDoseParameters(obj.getDoseParameters()/pln.numOfFractions);
        else
           obj = obj.setDoseParameters(obj.getDoseParameters());
        end
        cst{i,6}{j} = obj;        
    end
end

% resizing cst to dose cube resolution 
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,  dij.ctGrid.y,  dij.ctGrid.z,...
                                 dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

% find target indices and described dose(s) for weight vector
% initialization
V          = [];
doseTarget = [];
ixTarget   = [];

for i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}{1}];
        
        %Iterate through objectives/constraints
        fDoses = [];
        for fObjCell = cst{i,6}
            dParams = fObjCell{1}.getDoseParameters();
            %Don't care for Inf constraints
            dParams = dParams(isfinite(dParams));
            %Add do dose list
            fDoses = [fDoses dParams];
        end
                
        doseTarget = [doseTarget fDoses];
        ixTarget   = [ixTarget i*ones(1,length(fDoses))];
    end
end
[doseTarget,i] = max(doseTarget);
ixTarget       = ixTarget(i);

   
% calculate initial beam intensities wInit
matRad_cfg.dispInfo('Estimating initial weights... ');
if isfield(dij,'numOfModalities')
    numOfModalities = dij.numOfModalities;
else
    numOfModalities = 1;
end

bioModels = [pln.originalPlans.bioParam];

if exist('wInit','var')
    %do nothing as wInit was passed to the function
    matRad_cfg.dispInfo('chosen provided wInit!\n');   

else

    doseTarget = doseTarget/numOfModalities;
    wInit = [];
    % loop over all modalities
    for modality = 1: numOfModalities

        if isfield(dij,'original_Dijs')   
            dijt = dij.original_Dijs{modality};
            currentBioModel = pln.originalPlans(modality).bioParam;
        else 
            dijt = dij;
        end

        wOnes          = ones(dijt.totalNumOfBixels,1);
        wt             = zeros(dijt.totalNumOfBixels,1);
        
        if strcmp(currentBioModel.model,'constRBE') % If there is constRBE, the constRBE BP will be instantiated, and we need RBE field in dij.
            
            % check if a constant RBE is defined - if not use 1.1
            if ~isfield(dijt,'RBE')
                switch pln.originalPlans(modality).radiationMode
                    case 'protons'
                        dijt.RBE = 1.1;
                    case 'photons'
                        dijt.RBE = 1;
                end
            end

            doseTmp = dijt.physicalDose{1}*wOnes;
            bixelWeight =  (doseTarget)/(dijt.RBE * mean(doseTmp(V)));
            wt       = wOnes * bixelWeight;
            matRad_cfg.dispInfo('chosen uniform weight of %f!\n',bixelWeight);

        elseif any(strcmp(pln.bioParam.quantityOpt, {'RBExD', 'effect'})) %currentBioModel.bioOpt %%%% here we could have photons with RBExD that have bioOpt = 0;
            
            % retrieve photon LQM parameter
            [ax,bx] = matRad_getPhotonLQMParameters(cst,dijt.doseGrid.numOfVoxels,1);


            if ~isequal(dijt.ax(dijt.ax~=0),ax(dijt.ax~=0)) || ...
                    ~isequal(dijt.bx(dijt.bx~=0),bx(dijt.bx~=0))
                matRad_cfg.dispError('Inconsistent biological parameter - please recalculate dose influence matrix!\n');
            end


            for i = 1:size(cst,1)
                for j = 1:size(cst{i,6},2)
                    if strcmp(pln.radiationMode, 'MixMod')
                        DoseParameters = cst{i,6}{j}.getDoseParameters()./sum([dij.STfractions{:}]);
                    else
                        DoseParameters = cst{i,6}{j}.getDoseParameters();
                    end
                    % check if prescribed doses are in a valid domain
                    if any(DoseParameters > 5) && isequal(cst{i,3},'TARGET')
                        matRad_cfg.dispWarning('Reference dose > 10 Gy[RBE] for target. Biological optimization outside the valid domain of the base data. Reduce dose prescription or use more fractions.\n');
                    end

                end
            end

            dijt.ixDose  = dijt.bx~=0;
            dij.original_Dijs{modality}.ixDose = dijt.ixDose;

            if isequal(pln.bioParam.quantityOpt,'effect') %here check on the mixMod bioModel, the quantityOpt should be consistent

                effectTarget = cst{ixTarget,5}.alphaX * doseTarget + cst{ixTarget,5}.betaX * doseTarget^2;
                bixelNum = 1;
                SelectedBixels = [bixelNum:bixelNum+dijt.totalNumOfBixels*dij.numOfSTscen(modality)-1];
                wIdx = reshape(wOnes(SelectedBixels),[dijt.totalNumOfBixels,dij.numOfSTscen(modality)]);
                aTmp = dijt.mAlphaDose{1}*wIdx(:,1);
                bTmp = dijt.mSqrtBetaDose{1}*wIdx(:,1);

                p = sum(aTmp(V)) / sum(bTmp(V).^2);
                q = -(effectTarget * length(V)) / sum(bTmp(V,1).^2);

                wt(SelectedBixels)        = -(p/2) + sqrt((p^2)/4 - q) * wOnes(SelectedBixels);

                bixelNum = bixelNum + dijt.totalNumOfBixels;


            elseif isequal(pln.bioParam.quantityOpt,'RBExD')

                %pre-calculations
                dijt.gamma             = zeros(dijt.doseGrid.numOfVoxels,1);
                dijt.gamma(dijt.ixDose) = dijt.ax(dijt.ixDose)./(2*dijt.bx(dijt.ixDose));
                dij.original_Dijs{modality}.gamma = dijt.gamma;
                bixelNum = 1;

                SelectedBixels = [bixelNum:bixelNum+dijt.totalNumOfBixels*dij.numOfSTscen(modality)-1];
                wIdx = reshape(wOnes(SelectedBixels),[dijt.totalNumOfBixels,dij.numOfSTscen(modality)]);

                % calculate current effect in target
                aTmp = dijt.mAlphaDose{1}*wIdx(:,1);
                bTmp = dijt.mSqrtBetaDose{1} * wIdx(:,1);
                doseTmp = dijt.physicalDose{1}*wIdx(:,1);
                CurrEffectTarget = aTmp(V,1) + bTmp(V,1).^2;
                TolEstBio        = 1.2;
                % calculate maximal RBE in target
                maxCurrRBE = max(-cst{ixTarget,5}.alphaX + sqrt(cst{ixTarget,5}.alphaX^2 + ...
                    4*cst{ixTarget,5}.betaX.*CurrEffectTarget)./(2*cst{ixTarget,5}.betaX*doseTmp(V)));
                wt(SelectedBixels)    =  ((doseTarget)/(TolEstBio*maxCurrRBE*max(doseTmp(V))))* wOnes(SelectedBixels);
                bixelNum = bixelNum + dijt.totalNumOfBixels;

            
            end
            matRad_cfg.dispInfo('chosen weights adapted to biological dose calculation!\n');
        else
            bixelNum = 1;
            
            SelectedBixels = [bixelNum:bixelNum+dijt.totalNumOfBixels*dij.numOfSTscen(modality)-1];
            wIdx = reshape(wOnes(SelectedBixels),[dijt.totalNumOfBixels,dij.numOfSTscen(modality)]);

            doseTmp = dijt.physicalDose{1}*wIdx;
            bixelWeight =  (doseTarget)/mean(doseTmp(V));
            wt(SelectedBixels)       = wOnes(SelectedBixels) * bixelWeight;

            matRad_cfg.dispInfo('chosen uniform weight of %f!\n',bixelWeight);

        
        end
        
        wInit = [wInit; wt./max(wt)];            
    end
end


%% calculate probabilistic quantities for probabilistic optimization if at least
% one robust objective is defined 

%Check how to use 4D data
if isfield(pln,'propOpt') && isfield(pln.propOpt,'scen4D')
    scen4D = pln.propOpt.scen4D;
else
    scen4D = 1; %Use only first 4D scenario for optimization
end

%If "all" provided, use all scenarios
if isequal(scen4D,'all')

   scen4D = 1:size(dij.physicalDose,1);
end

if strcmp(pln.radiationMode, 'MixMod') && (scen4D ~= 1)
   scen4D = 1;
   matRad_dispWarning('For the time being, scen4D is disabled for MixMod optimization. \n');
end

linIxDIJ = [];
% loop over all modalities
for modality = 1: numOfModalities
    if isfield(dij,'original_Dijs')   
        dijt = dij.original_Dijs{modality};
    else 
        dijt = dij;
    end
     linIxDIJ = [linIxDIJ , find(~cellfun(@isempty,dijt.physicalDose(scen4D,:,:)))'];

     FLAG_CALC_PROB = false;
     FLAG_ROB_OPT   = false;

      for i = 1:size(cst,1)
          for j = 1:numel(cst{i,6})
             if strcmp(cst{i,6}{j}.robustness,'PROB') && numel(linIxDIJ{k}) > 1
                FLAG_CALC_PROB = true;
             end
             if ~strcmp(cst{i,6}{j}.robustness,'none') && numel(linIxDIJ{k}) > 1
                FLAG_ROB_OPT = true;
             end
          end
      end


   if FLAG_CALC_PROB
         if stcmp(pln.radiationMode, 'MixMod')

            [dij.original_Dijs{modality}] = matRad_calculateProbabilisticQuantities(dijt,cst,pln);
         else
            [dij] = matRad_calculateProbabilisticQuantities(dij,cst,pln);            
         end
   end
end

% set optimization options
if ~FLAG_ROB_OPT || FLAG_CALC_PROB     % if multiple robust objectives are defined for one structure then remove FLAG_CALC_PROB from the if clause
   ixForOpt = scen4D;
else
   ixForOpt = linIxDIJ;
end

switch pln.bioParam.quantityOpt
    case 'effect'
        backProjection = matRad_EffectProjection;
    case 'RBExD'
        %Capture special case of constant RBE
        protonPlnIdx = find(strcmp({pln.originalPlans.radiationMode}, 'protons'));

        if strcmp(pln.originalPlans(protonPlnIdx).bioParam.model, 'constRBE') %If proton has constRBE, photon does not but calcCombi dose computes the dij.RBE=1 anyway. TODO: implement this correctly
            backProjection = matRad_ConstantRBEProjection;
        else
            backProjection = matRad_VariableRBEProjection;
        end
    case 'physicalDose'
        backProjection = matRad_DoseProjection;
    otherwise
        warning(['Did not recognize bioloigcal setting ''' pln.probOpt.bioOptimization '''!\nUsing physical dose optimization!']);
        backProjection = matRad_DoseProjection;
end

%Give scenarios used for optimization
backProjection.scenarios    = ixForOpt;

backProjection.scenarioProb = [pln.multScen.scenProb];

if strcmp(pln.radiationMode, 'MixMod')
   optiProb = matRad_OptimizationProblemMixMod(backProjection);
   optiProb.nFractions = pln.numOfFractions;
else
   
   optiProb = matRad_OptimizationProblem(backProjection);
end
optiProb.quantityOpt = pln.bioParam.quantityOpt;
if isfield(pln,'propOpt') && isfield(pln.propOpt,'useLogSumExpForRobOpt')
    optiProb.useLogSumExpForRobOpt = pln.propOpt.useLogSumExpForRobOpt;
end

if dij.precon
    dij = matRad_mixModPreconditioner(dij);
end
%Get Bounds

if ~isfield(pln.propOpt,'boundMU')
    pln.propOpt.boundMU = 0;
end 
    % NOTE: for MixMod, as it is should not work now
%     if pln.propOpt(i).boundMU
%         if (isfield(dij,'minMU') || isfield(dij,'maxMU')) && ~isfield(dij,'numParticlesPerMU')
%             matRad_cfg.dispWarning('Requested MU bounds but number of particles per MU not set! Bounds will not be enforced and standard [0,Inf] will be used instead!');
%         elseif ~isfield(dij,'minMU') && ~isfield(dij,'maxMU')
%             matRad_cfg.dispWarning('Requested MU bounds but machine bounds not defined in dij.minMU & dij.maxMU! Bounds will not be enforced and standard [0,Inf] will be used instead!');
%         else
%             if isfield(dij,'minMU')
%                 optiProb.minimumW = dij.numParticlesPerMU .* dij.minMU / 1e6;
%                 matRad_cfg.dispInfo('Using lower MU bounds provided in dij!\n')
%             end
%     
%             if isfield(dij,'maxMU')
%                 optiProb.maximumW = dij.numParticlesPerMU .* dij.maxMU / 1e6;
%                 matRad_cfg.dispInfo('Using upper MU bounds provided in dij!\n')
%             end
%         end
%     else
%         matRad_cfg.dispInfo('Using standard MU bounds of [0,Inf]!\n')
%     end

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
       
optimizer = optimizer.optimize(wInit,optiProb,dij,cst);

wOpt = optimizer.wResult;
info = optimizer.resultInfo;
bxidx = 1;

for mod = 1: pln.numOfModalities

    wt = [];
    % split the w for current modality
    STrepmat = (~dij.spatioTemp(mod) + dij.spatioTemp(mod)*dij.numOfSTscen(mod));
    wt = reshape(wOpt(bxidx: bxidx+STrepmat*dij.original_Dijs{mod}.totalNumOfBixels-1),[dij.original_Dijs{mod}.totalNumOfBixels,STrepmat]);
    
    resultGUI{mod} = matRad_calcCubes(wt,dij.original_Dijs{mod});
    if isfield(dij,'preconW') && dij.precon
       wt = wt.*dij.preconW(mod);     
    end
    resultGUI{mod}.wUnsequenced = wt;
    resultGUI{mod}.usedOptimizer = optimizer;
    resultGUI{mod}.info = info;
    bxidx = bxidx + STrepmat*dij.original_Dijs{mod}.totalNumOfBixels;
end
%Robust quantities
if FLAG_ROB_OPT || numel(ixForOpt) > 1   
    Cnt = 1;
    for i = find(~cellfun(@isempty,dij.physicalDose))'
        tmpResultGUI = matRad_calcCubes(wOpt,dij,i);
        resultGUI.([pln.bioParam.quantityVis '_' num2str(Cnt,'%d')]) = tmpResultGUI.(pln.bioParam.quantityVis);
        Cnt = Cnt + 1;
    end
end

% unblock mex files
clear mex
