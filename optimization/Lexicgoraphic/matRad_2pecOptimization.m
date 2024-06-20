function [resultGUIs,resultGUIs2,cst1,cst2,PriorityList2]=  matRad_2pecOptimization(PriorityList1,dij,cst,pln,wInit)
    % Lexicographic optimization using the 2pec method as introducded by
    % Breedveld et al
    % 
    % call
    %
    % input
    %   PriorityList1:   Initial structure to store all objectives
    %   dij:            matRad dij struct
    %   cst:            modified matRad cst struct. 
    %   pln:            matRad pln struct
    %   wInit:          (optional) custom weights to initialize problems
    %
    % output
    %   resultGUIs:     Plans created in the first phase
    %   resultGUIs2:    Plans created in teh second phase
    %   cst1:           Modified cst after the first phase
    %   cst2:           Modified cst as after the second phase
    %   PriorityList2:  Final priority list for lexicographic optimization
    %
    % References
    %   https://iopscience.iop.org/article/10.1088/0031-9155/54/23/011
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024 the matRad development team. 
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

    PriorityList1.adaptToFractionSize(pln.numOfFractions);

    if exist('wInit','var') 
        [dij,cst,pln,wInit,optiProb] = matRad_initOptimization(dij,cst,pln,wInit);
    else
        [dij,cst,pln,wInit,optiProb] = matRad_initOptimization(dij,cst,pln);
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

    resultGUIs = {}; % resultGUIs for first phase
    resultGUIs2 = {}; % resultGUIs for second phase
    PriorityList2 = matRad_PriorityList2(); % create empty priority list
    PriorityList2.ConstraintList = PriorityList1.ConstraintList; %same constraints
    
    %initial update for cst with all objectives and constraints for overlap
    cst = PriorityList1.generateFullCst(cst);
    
    %generate the actual csts that get modified throughout the algorithm
    cst1 = PriorityList1.generateConstraintCst(cst);
    cst2 = cst1;

    matRad_cfg.propOpt.defaultAccChangeTol = 1e-6;
    %add first objective(s) to cst
    cst1 =  PriorityList1.modifyCst(cst1);
    while PriorityList1.numOfObj <= numel(PriorityList1.GoalList)

        optiProb.extractObjectivesAndConstraintsFromCst(cst1);
        

        %check if objective can be skipped
        [objectives,FastCalc] = PriorityList1.fastObjectiveCalc(dij,cst1,optiProb,wInit);

        if ~FastCalc % can the objective be skipped?
            %If not optimize the current objective
            optimizer = optimizer.optimize(wInit,optiProb,dij,cst1);
            
            wOpt = optimizer.wResult;
            info = optimizer.resultInfo;
            
            resultGUI = matRad_calcCubes(wOpt,dij);
            resultGUI.wUnsequenced = wOpt;
            resultGUI.usedOptimizer = optimizer;
            resultGUI.info = info;
            resultGUI.optiProb = optiProb;
            objectives = matRad_objectiveFunctions(optiProb,wOpt,dij,cst1);
    
            wInit = wOpt;

        else %objectives can be met
            resultGUI = matRad_calcCubes(wInit,dij);
            resultGUI.wUnsequenced = wInit;
            resultGUI.optiProb = optiProb;
            objectives = matRad_objectiveFunctions(optiProb,wInit,dij,cst1);
        end

        resultGUI.objectives = objectives;  %add objectives 

        resultGUIs{end+1} = resultGUI;
        % exchange previously optimized objectives to constraints and remove from priorityList
        [cst1,cst2,PriorityList2] = PriorityList1.updateStep(cst1,cst2,PriorityList2,objectives); 
        
        %add next objective
        if PriorityList1.numOfObj <= numel(PriorityList1.GoalList)
            cst1 = PriorityList1.modifyCst(cst1);
        end
    end

    %step2
    matRad_cfg.propOpt.defaultAccChangeTol = 1e-6;
    cst2 = PriorityList2.modifyCst(cst2); %should set objective in appropriate spot -> use VOIIdx
    while PriorityList2.numOfObj <= numel(PriorityList2.GoalList)
        optiProb.extractObjectivesAndConstraintsFromCst(cst2);
        optimizer = optimizer.optimize(wInit,optiProb,dij,cst2);
            
        wOpt = optimizer.wResult;
        info = optimizer.resultInfo;
        
        resultGUI = matRad_calcCubes(wOpt,dij);
        resultGUI.wUnsequenced = wOpt;
        resultGUI.usedOptimizer = optimizer;
        resultGUI.info = info;
        resultGUI.optiProb = optiProb;
        objectives = matRad_objectiveFunctions(optiProb,wOpt,dij,cst2);
        resultGUI.objectives = objectives;
        wInit = wOpt;
    

        resultGUIs2{end+1} = resultGUI;
        [cst2] = PriorityList2.updateStep(cst2,resultGUI.objectives); 

        if PriorityList2.numOfObj <= numel(PriorityList2.GoalList)
            cst2 = PriorityList2.modifyCst(cst2);
        end

    end
end