function targetDoseInfo = matRad_getTargetReferenceDoses(cst,pln)
% matRad_getTargetReferenceDoses resolves target reference doses from the cst
%
% call
%   targetDoseInfo = matRad_getTargetReferenceDoses(cst,pln)
%
% input
%   cst:                matRad cst cell array
%   pln:                matRad pln struct
%
% output
%   targetDoseInfo:     struct array with one entry per TARGET VOI. Each
%                       entry contains cstIndex, name, and refDose in dose
%                       per fraction.
%
% References
%   -
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

isTarget = cellfun(@(voiType) isequal(voiType,'TARGET'),cst(:,3));
targetRows = find(isTarget);
targetDoseInfo = repmat(struct('cstIndex',[],'name',[],'refDose',[]), ...
    numel(targetRows),1);

for refTarget = 1:numel(targetRows)
    cstIdx = targetRows(refTarget);

    refDose = getReferenceDoseFromObjectives(cst{cstIdx,6});
    if isfinite(refDose)
        refDose = refDose/pln.numOfFractions;
    else
        matRad_cfg.dispWarning('Target %s has no objective that penalizes underdosage. Reference dose unavailable.\n',cst{cstIdx,2});
    end

    targetDoseInfo(refTarget).cstIndex = cstIdx;
    targetDoseInfo(refTarget).name = cst{cstIdx,2};
    targetDoseInfo(refTarget).refDose = refDose;
end

end

function refDose = getReferenceDoseFromObjectives(objectives)
matRad_cfg = MatRad_Config.instance();
refDose = inf;

if isempty(objectives)
    return;
end

if isstruct(objectives)
    objectives = num2cell(arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,objectives));
end

for runObjective = 1:numel(objectives)
    obj = objectives{runObjective};
    if ~isa(obj,'matRad_DoseOptimizationFunction')
        try
            obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
        catch ME
            matRad_cfg.dispWarning('Objective/Constraint not valid!\n%s',ME.message);
            continue;
        end
    end

    isUnderdoseObjective = isa(obj,'DoseObjectives.matRad_SquaredDeviation') || ...
        isa(obj,'DoseObjectives.matRad_SquaredUnderdosing') || ...
        isa(obj,'DoseObjectives.matRad_MinDVH');

    if isUnderdoseObjective
        doseParameters = obj.getDoseParameters();
        doseParameters = doseParameters(isfinite(doseParameters));
        if ~isempty(doseParameters)
            refDose = min([refDose; doseParameters(:)]);
        end
    end
end
end
