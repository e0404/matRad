function assigned = matRad_recursiveFieldAssignment(assignTo,reference,fieldChangedWarningMessage,fieldname)
% matRad recursive field assignment tool
%   This function recursively assigns fields from one structure to another. If both 'assignTo' and 'reference' are structures,
% it will recurse into their fields. If a field in 'assignTo' is a structure and its corresponding field in 'reference' is not,
% or vice versa, a warning message is displayed. The function also handles the case where 'assignTo' or 'reference' are not structures,
% directly assigning the values. Custom warning messages can be specified for overwriting fields.
%
% call
%   assigned = matRad_recursiveFieldAssignment(assignTo,reference,fieldChangedWarningMessage,fieldname)
%
% input
%   assignTo:                      The initial structure to which the fields are to be assigned.
%   reference:                     The structure containing the fields and values to be assigned to 'assignTo'.
%   fieldChangedWarningMessage:    Optional. A message to display if a field is overwritten. If not provided, no message is displayed.
%   fieldname:                     Optional. The name of the current field being processed. Used for generating specific warning messages.
%
% output
%   assigned:                      The structure 'assignTo' after assigning the fields from 'reference'.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3 || isempty(fieldChangedWarningMessage)
    fieldChangedWarningMessage = '';
else
    if ~isempty(fieldChangedWarningMessage) && fieldChangedWarningMessage(end) ~= ':'
        fieldChangedWarningMessage = [fieldChangedWarningMessage ':'];
    end
end

if nargin < 4
    fieldname = '';
end

% If the reference is a struct, we need to recurse into it
if isstruct(assignTo) && isstruct(reference)
    %First make sure the output has all the fields of the assignTo structure
    assigned = assignTo;
    
    %Now iterate over all fields provided in the reference structure
    fields = fieldnames(reference);
    for i = 1:numel(fields)
        field = fields{i};
        if ~isfield(assignTo,field)
            assigned.(field) = reference.(field);
        else
            assigned.(field) = matRad_recursiveFieldAssignment(assignTo.(field),reference.(field),fieldChangedWarningMessage,field);
        end
    end
% If the reference is not a struct, we can assign it directly. 
% However, we need to check if the assignTo is a struct and if we need to warn the user of overwriting a struct with a non-struct
elseif isstruct(assignTo) && ~isstruct(reference)        
    if ~isempty(fieldname)
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispWarning([fieldChangedWarningMessage 'Field ''%s'' is a struct but will be overwritten by a ''%s!'''],fieldname,class(reference));
    end
    assigned = reference;
% If the assignTo is not a struct, we can assign it directly.
% However, we need to check if the reference is a struct and if we need to warn the user of overwriting a non-struct with a struct
elseif ~isstruct(assignTo) && isstruct(reference)
    if ~isempty(fieldname)
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispWarning([fieldChangedWarningMessage 'Field ''%s'' is not a struct but will be overwritten by a struct!'],fieldname);
    end
    assigned = reference;
else
    if ~isempty(fieldname) && ~isempty(fieldChangedWarningMessage)
        matRad_cfg = MatRad_Config.instance();
        if isequal(class(assignTo),class(reference))
            matRad_cfg.dispWarning([fieldChangedWarningMessage 'Field ''%s'' will be overwritten by reference value!'],fieldname);
        else
            matRad_cfg.dispWarning([fieldChangedWarningMessage 'Field ''%s'' is supposed to be a %s but will be overwritten by a ''%s!'''],fieldname,class(assignTo),class(reference));
        end
    end
    assigned = reference;
end

end