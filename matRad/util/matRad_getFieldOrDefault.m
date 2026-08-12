function val = matRad_getFieldOrDefault(s, fieldName, defaultVal, warningMessage, deprecation)
% matRad_getFieldOrDefault Reads s.(fieldName), falling back to a default
%   if the struct or the field doesn't exist, optionally displaying a
%   warning (or deprecation warning) when the fallback is used.
%
% call:
%   val = matRad_getFieldOrDefault(s, fieldName, defaultVal)
%   val = matRad_getFieldOrDefault(s, fieldName, defaultVal, warningMessage)
%   val = matRad_getFieldOrDefault(s, fieldName, defaultVal, warningMessage, deprecation)
%
% input:
%   s:               struct to read the field from (may be missing/empty/non-struct)
%   fieldName:       name of the field to read
%   defaultVal:      value returned if the field is not present
%   warningMessage:  optional. If given (non-empty), displayed via
%                    matRad_cfg.dispWarning (or dispDeprecationWarning, see
%                    below) whenever the field is missing and defaultVal is
%                    used instead
%   deprecation:     optional. If true, warningMessage is displayed via
%                    matRad_cfg.dispDeprecationWarning instead of
%                    dispWarning. Default: false
%
% output:
%   val:             s.(fieldName) if present, otherwise defaultVal
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

if nargin < 4
    warningMessage = '';
end
if nargin < 5 || isempty(deprecation)
    deprecation = false;
end

if isstruct(s) && isfield(s, fieldName)
    val = s.(fieldName);
else
    val = defaultVal;
    if ~isempty(warningMessage)
        matRad_cfg = MatRad_Config.instance();
        if deprecation
            matRad_cfg.dispDeprecationWarning(warningMessage);
        else
            matRad_cfg.dispWarning(warningMessage);
        end
    end
end

end
