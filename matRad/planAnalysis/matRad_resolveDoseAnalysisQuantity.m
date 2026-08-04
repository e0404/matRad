function quantity = matRad_resolveDoseAnalysisQuantity(resultGUI, pln, requestedQuantity)
% matRad_resolveDoseAnalysisQuantity resolves the resultGUI dose field.
%
% call
%   quantity = matRad_resolveDoseAnalysisQuantity(resultGUI,pln)
%   quantity = matRad_resolveDoseAnalysisQuantity(resultGUI,pln,requestedQuantity)
%
% input
%   resultGUI:           matRad resultGUI struct
%   pln:                 matRad pln struct
%   requestedQuantity:   optional explicit resultGUI dose field
%
% output
%   quantity:            selected resultGUI dose field
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

if nargin < 3
    requestedQuantity = '';
end

quantity = matRad_charIfString(requestedQuantity);
if ~isempty(quantity)
    matRad_validateQuantity(resultGUI, quantity);
    return
end

quantity = matRad_getPlanQuantity(pln, 'quantityOpt');
if ~isempty(quantity) && isfield(resultGUI, quantity)
    return
end

if isfield(resultGUI, 'analysisQuantity') && ...
        isfield(resultGUI, matRad_charIfString(resultGUI.analysisQuantity))
    quantity = matRad_charIfString(resultGUI.analysisQuantity);
elseif isfield(resultGUI, 'RBExDose')
    quantity = 'RBExDose';
else
    quantity = 'physicalDose';
end

matRad_validateQuantity(resultGUI, quantity);

end

function quantity = matRad_getPlanQuantity(pln, fieldName)
quantity = '';
if isfield(pln, 'propOpt') && isfield(pln.propOpt, fieldName) && ...
        ~isempty(pln.propOpt.(fieldName))
    quantity = matRad_charIfString(pln.propOpt.(fieldName));
end
end

function value = matRad_charIfString(value)
if isstring(value)
    value = char(value);
end
end

function matRad_validateQuantity(resultGUI, quantity)
if ~isfield(resultGUI, quantity)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Unknown quantity ''%s'' to analyse!', quantity);
end
end
