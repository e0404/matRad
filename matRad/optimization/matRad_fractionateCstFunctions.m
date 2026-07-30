function cst = matRad_fractionateCstFunctions(cst, numOfFractions)
% matRad function to prepare the optimization functions stored in the cst
% for a fractionated optimization
%
% Promotes every objective/constraint in cst{:,6} to an optimization function
% object (accepting the legacy struct formats) and rescales the dose
% parameters of the dose related ones from the prescribed total dose to the
% dose per fraction. Functions that do not act on a dose related quantity
% (e.g. fluence objectives) are promoted but left untouched.
%
% call:
%   cst = matRad_fractionateCstFunctions(cst,numOfFractions)
%
% input:
%   cst:              matRad cst struct
%   numOfFractions:   number of fractions the plan is delivered in
%
% output:
%   cst:              matRad cst struct with objects in cst{:,6} whose dose
%                     parameters refer to a single fraction
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

matRad_cfg = MatRad_Config.instance();

for i = 1:size(cst, 1)
    % compatibility layer for the old objective format, where all functions
    % of a structure were stored as a struct array instead of a cell array
    if isstruct(cst{i, 6})
        cst{i, 6} = arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct, cst{i, 6}, 'UniformOutput', false);
    end

    for j = 1:numel(cst{i, 6})

        obj = cst{i, 6}{j};

        % In case it is a default saved struct, convert to object
        % Also intrinsically checks that we have a valid optimization
        % objective or constraint function in the end
        if ~isa(obj, 'matRad_OptimizationFunction')
            try
                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
            catch
                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!', i, j);
            end
        end

        % only dose related functions carry dose parameters to fractionate
        if isa(obj, 'matRad_DoseOptimizationFunction')
            obj = obj.setDoseParameters(obj.getDoseParameters() / numOfFractions);
        end

        cst{i, 6}{j} = obj;
    end
end

end
