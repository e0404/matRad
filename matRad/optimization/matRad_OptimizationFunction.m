classdef (Abstract) matRad_OptimizationFunction
    % matRad_OptimizationFunction. Superclass for all optimization functions.
    % This is the common ancestor of every objective and constraint, no matter
    % which quantity it acts on. It only provides what is genuinely independent
    % of that quantity: the parameter description used by the GUI, struct
    % serialization and instantiation from a struct.
    %
    % The quantity specific bases derive from it::
    %
    %   matRad_DoseOptimizationFunction     functions of the dose (adds the dose
    %                                       parameter accessors and robustness)
    %   matRad_FluenceOptimizationFunction  functions of the fluence vector
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

    properties (Abstract, Constant)
        name                % Display name of the function. Needs to be implemented in sub-classes.
        parameterNames      % Cell array of display names of the parameters. Needs to be implemented in sub-classes.
        parameterTypes      % Cell array of parameter types. Valid types are 'dose', 'numeric', or a cell list of string options.
    end

    properties (Abstract, Access = public)
        parameters
    end

    methods

        function obj = matRad_OptimizationFunction(dataStruct)
            if nargin > 0 && ~isempty(dataStruct) && isstruct(dataStruct)
                obj = assignCommonPropertiesFromStruct(obj, dataStruct);
            end
        end

        % Overload the struct function to return a struct with the general
        % information of the objective / constraint
        function s = struct(obj)
            s.className = class(obj);
            s.parameters = obj.parameters;
        end

    end

    methods (Access = private)

        function obj = assignCommonPropertiesFromStruct(obj, s)
            for fn = fieldnames(s)'    % enumerate fields
                try
                    obj.(fn{1}) = s.(fn{1});   % and copy
                catch
                    continue
                    % Do Nothing here
                    % warning('Could not copy field %s', fn{1});
                end
            end
        end

    end

    methods (Static)

        % creates an optimization function from a struct
        function obj = createInstanceFromStruct(s)
            try
                % Create objective / constraint from class name
                obj = eval([s.className '(s)']);

                env = matRad_getEnvironment();

                % Workaround for Octave which has a problem assigning
                % properties in superclass
                if strcmp(env, 'OCTAVE')
                    obj = assignCommonPropertiesFromStruct(obj, s);
                end

            catch ME
                error(['Could not instantiate Optimization Function: ' ME.message]);
            end
        end

    end
end
