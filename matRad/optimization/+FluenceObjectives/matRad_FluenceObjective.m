classdef (Abstract) matRad_FluenceObjective < matRad_FluenceOptimizationFunction
    % matRad_FluenceObjective: Interface for optimization objectives acting on
    %   the fluence
    %   This abstract base class provides the structure of objectives that are a
    %   function of the bixel weight vector w itself, not of the dose. They take
    %   no part in the dose backprojection / chain rule - the gradient they
    %   return is already a gradient w.r.t. the weights and is added to the
    %   weight gradient directly. Implementations can be found in the
    %   FluenceObjectives package.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Abstract, Access = public)
        penalty                 % Optimization penalty, applied by the objective wrapper
    end

    methods (Access = public)

        function obj = matRad_FluenceObjective(varargin)
            % default initialization from struct (parameters & penalty)
            obj@matRad_FluenceOptimizationFunction(varargin{:});
        end

        % Overloads the struct function to add objective related information
        function s = struct(obj)
            s = struct@matRad_FluenceOptimizationFunction(obj);
            s.penalty = obj.penalty;
        end

    end

    % These should be abstract methods, however Octave can't parse them. As
    % soon as Octave is able to do this, they should be made abstract again
    methods % (Abstract)

        % returns the objective function value for the given fluence vector
        function fFluence = computeFluenceObjectiveFunction(obj, w) %#ok<INUSD>
            error('Function needs to be implemented!');
        end

        % returns the fluence gradient for the given fluence vector
        function fFluenceGrad = computeFluenceObjectiveGradient(obj, w) %#ok<INUSD>
            error('Function needs to be implemented!');
        end

    end

end
