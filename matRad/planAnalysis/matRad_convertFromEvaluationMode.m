function [baseValue,evaluationMode,evaluationScale] = matRad_convertFromEvaluationMode(value,pln,evaluationMode)
% matRad_convertFromEvaluationMode Convert evaluation values to per-fraction.
%
% call
%   baseValue = matRad_convertFromEvaluationMode(value,pln,evaluationMode)
%   [baseValue,evaluationMode,evaluationScale] = matRad_convertFromEvaluationMode(value,pln,evaluationMode)
%
% input
%   value:             numeric value, array, or empty value expressed in
%                      evaluationMode
%   pln:               matRad pln struct
%   evaluationMode:    'perFraction' or 'total'
%
% output
%   baseValue:         value converted to the per-fraction evaluation base
%   evaluationMode:    normalized evaluation mode string
%   evaluationScale:   scalar conversion factor from base to evaluationMode
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

[~,evaluationMode,evaluationScale] = matRad_convertToEvaluationMode( ...
    1,pln,evaluationMode);
baseValue = value ./ evaluationScale;

end
