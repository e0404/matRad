function apertureInfo = matRad_upgradeApertureInfo(apertureInfo)
% matRad_upgradeApertureInfo Migrates a legacy apertureInfo struct to the
%   current field naming, so that aperture info saved by an older matRad
%   version (e.g. loaded from a .mat alongside a resultGUI) keeps working.
%
%   Renamed fields are moved to their current name and the legacy field is
%   removed, so downstream code only ever sees the current names. A
%   deprecation warning is issued once per renamed field. Structs that are
%   already up to date are returned unchanged and silently.
%
% call:
%   apertureInfo = matRad_upgradeApertureInfo(apertureInfo)
%
% input:
%   apertureInfo:   aperture weight and shape info struct, possibly using
%                   legacy field names
%
% output:
%   apertureInfo:   the same struct using the current field names
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

if ~isstruct(apertureInfo) || ~isfield(apertureInfo, 'beam')
    return
end

% Legacy -> current field names, per level of the struct. Only fields that
% were part of a released matRad version need an entry here: names that only
% ever existed on an unmerged branch cannot appear in anyone's saved data.
beamRenames = { ...
               'lim_l', 'limLeft'; ...
               'lim_r', 'limRight'};

apertureInfo.beam = matRad_upgradeFields(apertureInfo.beam, beamRenames, 'apertureInfo.beam');

end

function s = matRad_upgradeFields(s, renames, contextName)
% Applies a legacy -> current rename table to every element of the struct
% array s, warning once per renamed field rather than once per element.

if isempty(s)
    return
end

matRad_cfg = MatRad_Config.instance();

for r = 1:size(renames, 1)
    oldName = renames{r, 1};
    newName = renames{r, 2};

    if ~isfield(s, oldName)
        continue
    end

    matRad_cfg.dispDeprecationWarning('%s.%s is deprecated, using %s.%s instead.', ...
                                      contextName, oldName, contextName, newName);

    for i = 1:numel(s)
        % An explicit new field wins over the legacy one, so that a struct
        % carrying both (e.g. partially migrated by hand) is not corrupted.
        if ~isfield(s, newName) || isempty(s(i).(newName))
            s(i).(newName) = s(i).(oldName);
        end
    end

    s = rmfield(s, oldName);
end

end
