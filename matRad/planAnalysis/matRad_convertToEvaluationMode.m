function [convertedValue,evaluationMode,evaluationScale] = matRad_convertToEvaluationMode(rawValue,pln,evaluationMode)
% matRad_convertToEvaluationMode Convert per-fraction values for evaluation.
%
% call
%   convertedValue = matRad_convertToEvaluationMode(rawValue,pln,evaluationMode)
%   [convertedValue,evaluationMode,evaluationScale] = matRad_convertToEvaluationMode(rawValue,pln,evaluationMode)
%
% input
%   rawValue:          numeric value, array, or empty per-fraction value
%   pln:               matRad pln struct
%   evaluationMode:    'perFraction' or 'total'
%
% output
%   convertedValue:    rawValue converted for the requested evaluation mode
%   evaluationMode:    normalized evaluation mode string
%   evaluationScale:   scalar conversion factor applied to rawValue
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

if nargin < 3 || isempty(evaluationMode)
    evaluationMode = 'perFraction';
end

if isstring(evaluationMode)
    evaluationMode = char(evaluationMode);
end

switch lower(evaluationMode)
    case 'perfraction'
        evaluationMode = 'perFraction';
        evaluationScale = 1;
    case 'total'
        evaluationMode = 'total';
        evaluationScale = getNumOfFractions(pln);
    otherwise
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispError('Unknown evaluationMode "%s". Use "perFraction" or "total".',evaluationMode);
end

convertedValue = rawValue .* evaluationScale;

end

function numOfFractions = getNumOfFractions(pln)

numOfFractions = 1;
if isfield(pln,'numOfFractions') && ~isempty(pln.numOfFractions)
    numOfFractions = pln.numOfFractions;
end

if ~(isnumeric(numOfFractions) && isscalar(numOfFractions) && ...
        isfinite(numOfFractions) && numOfFractions > 0)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('pln.numOfFractions must be a positive finite scalar.');
end

end
