function [ix, reason] = matRad_resolveClassSelection(requested, availableClasses, allClasses, pln)
% matRad function to select a class (e.g. a dose engine) by its shortName and
% to explain why a requested class cannot be used
%
% The requested name is matched case-insensitively against the shortNames of
% the classes available for the plan. A cell array of names selects the first
% name that is available, in the order given, like the defaults in
% MatRad_Config. If no name matches, reason explains for every name whether it
% is unknown, does not support the radiation mode of the plan (according to its
% constant property possibleRadiationModes) or why it is not available, using
% the message returned by the isAvailable method of the class (#943).
%
% call:
%   [ix, reason] = matRad_resolveClassSelection(requested, availableClasses, allClasses, pln)
%
% input:
%   requested:          requested shortName as char or string, or a cell array
%                       of shortNames
%   availableClasses:   struct array of the classes available for the plan with
%                       the fields shortName and className (see
%                       matRad_identifyClassesByConstantProperties)
%   allClasses:         struct array of all classes (available or not) in the
%                       same format, or a function handle returning it. It is
%                       only evaluated if no requested name is available.
%   pln:                matRad pln struct passed to the isAvailable methods
%
% output:
%   ix:                 index into availableClasses of the selected class,
%                       empty if no requested name is available
%   reason:             explanation why no requested class can be used, empty
%                       if a class was selected
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

ix = [];
reason = '';

if ischar(requested) || isstring(requested) || iscellstr(requested)
    requestedNames = cellstr(requested);
else
    reason = 'the selection is neither a name nor a cell array of names';
    return
end

if isempty(availableClasses)
    availableNames = {};
else
    availableNames = {availableClasses(:).shortName};
end

for i = 1:numel(requestedNames)
    match = find(strcmpi(availableNames, requestedNames{i}), 1);
    if ~isempty(match)
        ix = match;
        return
    end
end

if isa(allClasses, 'function_handle')
    allClasses = allClasses();
end

reasons = cellfun(@(name) matRad_explainUnavailableClass(name, availableNames, allClasses, pln), requestedNames, 'UniformOutput', false);
reason = strjoin(reasons, '; ');

end

function reason = matRad_explainUnavailableClass(name, availableNames, allClasses, pln)
% Explains why the class with the given shortName is not available

match = find(strcmpi({allClasses(:).shortName}, name), 1);

if isempty(match)
    reason = sprintf('''%s'' is unknown (available: %s)', name, strjoin(availableNames, ', '));
    return
end

className = allClasses(match).className;

% The radiation mode is checked first, because the availability checks of
% some classes report misleading messages for other radiation modes
radiationModes = matRad_getPossibleRadiationModes(className);
if ~isempty(radiationModes) && ~any(strcmp(radiationModes, pln.radiationMode))
    reason = sprintf('''%s'' does not support radiation mode ''%s'' (supported: %s)', name, pln.radiationMode, strjoin(radiationModes, ', '));
    return
end

[available, msg] = matRad_callIsAvailable(className, pln);

if available
    % e.g. if the availability check of the class fails within the plan's
    % class list but not here
    reason = sprintf('''%s'' is not available for this plan', name);
elseif ~isempty(msg)
    reason = sprintf('''%s'' is not available: %s', name, msg);
elseif isfield(pln, 'machine')
    reason = sprintf('''%s'' does not support radiation mode ''%s'' with machine ''%s''', name, pln.radiationMode, pln.machine);
else
    reason = sprintf('''%s'' does not support radiation mode ''%s''', name, pln.radiationMode);
end

end

function [available, msg] = matRad_callIsAvailable(className, pln)
% Calls the static isAvailable method of a class, which may also return only
% the availability without a message. eval is used because str2func does not
% work on static methods in older Octave versions.

msg = '';
try
    [available, msg] = eval([className '.isAvailable(pln)']);
catch
    try
        available = eval([className '.isAvailable(pln)']);
    catch ME
        available = false;
        msg = sprintf('its availability check failed (%s)', ME.message);
    end
end

end

function radiationModes = matRad_getPossibleRadiationModes(className)
% Returns the constant property possibleRadiationModes of a class, or an empty
% cell if the class has none

try
    radiationModes = cellstr(eval([className '.possibleRadiationModes']));
catch
    radiationModes = {};
end

end
