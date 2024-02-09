classdef(Abstract) matRad_PriorityClass < handle
% matRad_PriorityClass creates the super class for PriorityLists used in
% the lexicographic approach
%
% References
%   -
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
    properties
        Priorities;
        ConstraintList;
        GoalList;
        slackVariable = 1.03;
        numOfObj;
    end

        
    methods
        function obj = matRad_PriorityClass()
            %add optional slack Variable adjustment
            obj.Priorities = []; % list storing all priorities (e.g [1,2,2,3])
            obj.GoalList = {};
            obj.ConstraintList = {};
            obj.numOfObj = 1;

        end
        
        function addObjective(obj,Priority,objective,goal,cstIdx)
        % function to add an objective to a priorityList
        %
        % input
        %   Priority:   priority of the objective in the priority list
        %   objective:  corresponding objective (either struct or class based)
        %   goal:       goal value associated with the objetive
        %   cstIdx:     Index of the VOI in the cst that this objective
        %               belongs to

            obj.Priorities = [obj.Priorities,Priority]; %add priority

            if ~isa(objective,'matRad_DoseOptimizationFunction')
                objective = matRad_DoseOptimizationFunction.createInstanceFromStruct(objective);          
            end
  
            obj.GoalList{end+1} = matRad_PriorityListObjective(objective,goal,cstIdx);
            %sort by priority
            [obj.Priorities,I] = sort(obj.Priorities);
            obj.GoalList = obj.GoalList(I);
        end

        function addConstraint(obj,constraint,cstIdx)
        % function to add an constraint to a priorityList
        %
        % input
        %   constraint: Constraint to be added (either struct or class based)
        %   cstIdx:     Index of the VOI in the cst that this objective
        %               belongs to
        %
            if ~isa(constraint,'matRad_DoseOptimizationFunction')
                constraint = matRad_DoseOptimizationFunction.createInstanceFromStruct(constraint);
            end

            obj.ConstraintList{end+1} = matRad_PriorityListConstraint(constraint,cstIdx);
        end

        function adaptToFractionSize(obj,numOfFractions)
        % function to adapt all constraints and objectives to the fraction
        % size level
            %adapt constraints to fraction size
            for i = 1:numel(obj.ConstraintList)
                constraint = obj.ConstraintList{i}.constraint;
                obj.ConstraintList{i}.constraint = constraint.setDoseParameters(constraint.getDoseParameters()/numOfFractions);
            end
            %adapt objectives dose and goal values to fraction size
            for i= 1:numel(obj.GoalList)
                objective = obj.GoalList{i}.objective;
                obj.GoalList{i}.objective = objective.setDoseParameters(objective.getDoseParameters()/numOfFractions);
                obj.GoalList{i}.goalValue = objective.adaptGoalToFraction(obj.GoalList{i}.goalValue,numOfFractions);
            end
        end
    %}

        function cst = generateConstraintCst(obj,cst)
            %generate a bare bone cst struct for prioritized optimization containing only the 
            % hard constraints
            cst(:,6) = cell(1);
            for i = 1:numel(obj.ConstraintList)
                cst{obj.ConstraintList{i}.cstIdx,6}{end+1} = obj.ConstraintList{i}.constraint;
            end
        end

        function [Priority,LexiObjectives] = nextPriority(obj) %get next element(s) in prioList
            Priority = obj.Priorities(obj.numOfObj); %get next highest priority
            LexiObjectives = {}; %cell that is filled with the objectives of next highest prio
            i = obj.numOfObj;
            % get all objectives for this Priority
            while i <= numel(obj.Priorities) && obj.Priorities(i) == Priority 
                LexiObjectives{end+1} = obj.GoalList{i};
                i = i + 1;
            end
        end

        function printPriorityList(obj,cst)
            %%%
            % Function to print out the current PriorityList in table form
            %%%
            %create data
            goals = [];
            VOINames= {};
            Objectives = {};
            AchievedValues = [];
            AchievedValues2 = [];
            for i = 1:numel(obj.GoalList)
                idx = obj.GoalList{i}.cstIdx;
                goals = [goals;obj.GoalList{i}.goalValue];
                VOINames = [VOINames;cst{idx,2}];
                Objectives = [Objectives;obj.GoalList{i}.objective.name];
                AchievedValues = [AchievedValues;obj.GoalList{i}.achievedValue];
                AchievedValues2 = [AchievedValues2;obj.GoalList{i}.achievedValue2];
            end
            Priority = obj.Priorities';
            T = table(Priority,VOINames,Objectives,goals,AchievedValues,AchievedValues2);
            fig = uifigure();
            uitable(fig,"Data",T);
        end
            

    end

    methods (Abstract)
        
        modifyCst
        updateStep  
    end

end