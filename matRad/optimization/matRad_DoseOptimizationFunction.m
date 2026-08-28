classdef (Abstract) matRad_DoseOptimizationFunction < matRad_OptimizationFunction
    % matRad_DoseOptimizationFunction. Superclass for dose objectives and
    % constraints.
    % This is the superclass for all objectives and constraints acting on the
    % dose, to enable easy one-line identification. It adds the dose parameter
    % accessors and the robustness setting to the generic
    % matRad_OptimizationFunction.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019-2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        robustness = 'none'     % Robustness setting -> may be removed from the DoseObjective class in a future release
    end

    methods

        function obj = matRad_DoseOptimizationFunction(varargin)
            obj@matRad_OptimizationFunction(varargin{:});
        end

        % Overload the struct function to return a struct with general
        % the objective / constraint
        function s = struct(obj)
            s = struct@matRad_OptimizationFunction(obj);
            s.robustness = obj.robustness;
        end

        function obj = set.robustness(obj, robustness)
            matRad_cfg = MatRad_Config.instance();
            if ischar(robustness) && any(strcmp(robustness, obj.availableRobustness()))
                obj.robustness = robustness;
            else
                matRad_cfg.dispError('Invalid robustness setting!');
            end
        end

    end

    % Helper methods
    methods (Access = public)

        function doseParams = getDoseParameters(obj)
            % get only the dose related parameters.
            ix = cellfun(@(c) isequal('dose', c), obj.parameterTypes);
            doseParams = [obj.parameters{ix}];
        end

        function obj = setDoseParameters(obj, doseParams)
            % set only the dose related parameters.
            ix = cellfun(@(c) isequal('dose', c), obj.parameterTypes);
            obj.parameters(ix) = num2cell(doseParams);
        end

    end

    methods (Static)

        % creates an optimization function from a struct
        function obj = createInstanceFromStruct(s)
            try
                % Check for old version of cst objectives / constraints and
                % convert if necessary
                if isfield(s, 'type')
                    s = matRad_DoseOptimizationFunction.convertOldOptimizationStruct(s);
                end
            catch ME
                error(['Could not instantiate Optimization Function: ' ME.message]);
            end

            % The generic instantiation lives in the common ancestor
            obj = matRad_OptimizationFunction.createInstanceFromStruct(s);
        end

        function rob = availableRobustness()
            rob = {'none'}; % By default, no robustness is available
        end

        function s = convertOldOptimizationStruct(oldStruct)
            % Converts old version objectives to current format
            switch oldStruct.type
                % Objectives
                case 'square deviation'
                    s.className = 'DoseObjectives.matRad_SquaredDeviation';
                    s.penalty = oldStruct.penalty;
                    s.parameters{1} = oldStruct.dose;

                case 'square overdosing'
                    s.className = 'DoseObjectives.matRad_SquaredOverdosing';
                    s.penalty = oldStruct.penalty;
                    s.parameters{1} = oldStruct.dose;

                case 'square underdosing'
                    s.className = 'DoseObjectives.matRad_SquaredUnderdosing';
                    s.penalty = oldStruct.penalty;
                    s.parameters{1} = oldStruct.dose;

                case 'min DVH objective'
                    s.className = 'DoseObjectives.matRad_MinDVH';
                    s.parameters{1} = oldStruct.dose;
                    s.parameters{2} = oldStruct.volume;

                case 'max DVH objective'
                    s.className = 'DoseObjectives.matRad_MaxDVH';
                    s.parameters{1} = oldStruct.dose;
                    s.parameters{2} = oldStruct.volume;

                case 'mean'
                    s.className = 'DoseObjectives.matRad_MeanDose';
                    s.parameters{1} = oldStruct.dose;

                case 'EUD'
                    s.className = 'DoseObjectives.matRad_EUD';
                    s.parameters{1} = oldStruct.dose;
                    s.parameters{2} = oldStruct.EUD;

                    % Constraints
                case  'max dose constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxDose';
                    s.parameters{1} = 0;
                    s.parameters{2} = oldStruct.dose;

                case  'min dose constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxDose';
                    s.parameters{1} = oldStruct.dose;
                    s.parameters{2} = Inf;

                case  'min mean dose constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxMeanDose';
                    s.parameters{1} = oldStruct.dose;
                    s.parameters{2} = Inf;

                case  'max mean dose constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxMeanDose';
                    s.parameters{1} = 0;
                    s.parameters{2} = oldStruct.dose;

                case  'min EUD constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxEUD';
                    s.parameters{1} = oldStruct.EUD;
                    s.parameters{2} = oldStruct.dose;
                    s.parameters{3} = Inf;

                case  'max EUD constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxEUD';
                    S.parameters{1} = oldStruct.EUD;
                    s.parameters{2} = 0;
                    s.parameters{3} = oldStruct.dose;

                case  'max DVH constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxDVH';
                    s.parameters{1} = oldStruct.dose;
                    s.parameters{2} = 0;
                    s.parameters{3} = oldStruct.volume;

                case  'min DVH constraint'
                    s.className = 'DoseConstraints.matRad_MinMaxDVH';
                    s.parameters{1} = oldStruct.dose;
                    s.parameters{2} = oldStruct.volume;
                    s.parameters{3} = 100;
                otherwise
                    ME = MException('optimization:ObjectCreationFailed', ...
                                    'Old versioned input struct / parameter invalid for creation of optimization function!');
                    throw(ME);
            end
        end

    end
end
