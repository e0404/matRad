function fileExists = matRad_checkMexFileExists(filename,linkOctave)
% Checks if a matching mex file exists, and can create a link
% if a matching custom, system specific precompiled octave mex file is found
%
% call:
%   matRad_checkMexFileExists(filename)
%   matRad_checkMexFileExists(filename,linkOctave)
%
% input:
%   filename:   name of the mex file (without extension)
%   linkOctave: (optional: default true) If set to true, the function will check
%               for a custom build mex file for octave with our custom extension
%               mexoct<version><system>. From Octave 10 on, <version> is the
%               major version only (e.g. mexoct10w64) because mex files link
%               against the stable liboctmex library and are compatible across
%               major versions; the newest file at or below the running major
%               version is used. Before Octave 10 the exact full version is used
%               (e.g. mexoct920w64). If such a file exists, it will create a link
%               to the file as filename.mex since octave does not use system-
%               specific extensions by default. If false, the function will only
%               check for an existing .mex file for octave.
% output:
%   fileExists: true if the mex file exists (or can be linked), and false
%               otherwise
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

if nargin < 2
    linkOctave = true;
end

%Explicitly check for matching mex file (id 3)
fileExists = (exist([filename '.' mexext],'file') == 3);

%For octave we ship experimental system-specific precompiled mex files
[env,ver] = matRad_getEnvironment();

if ~fileExists && matRad_cfg.isOctave && linkOctave

    %Check Architecture
    [~,maxArraySize] = computer();

    if maxArraySize > 2^31
        bitExt = '64';
    else
        bitExt = '32';
    end

    %Operating-system tag used in our custom extension
    if ispc
        osExt = 'w';
    elseif ismac
        osExt = 'mac';
    elseif isunix
        osExt = 'a';
    else
        %No file for unknown operating system
        fileExists = false;
        return;
    end

    %Determine which version tags to look for, newest first.
    %From Octave 10 on, mex files link against the stable liboctmex library
    %(instead of liboctave/liboctinterp), so they stay compatible across major
    %versions. We therefore ship one file per MAJOR version (mexoct10, mexoct11,
    %...) and let a given Octave pick the newest compatible file it can find by
    %iterating from its own major version down to 10. Before version 10 we keep
    %the exact full-version tag (e.g. mexoct920), as those files are not
    %cross-version compatible.
    majorVer = str2double(strtok(ver,'.'));
    if ~isnan(majorVer) && majorVer >= 10
        versionTags = arrayfun(@num2str, majorVer:-1:10, 'UniformOutput', false);
    else
        versionTags = {erase(ver,'.')};
    end

    %Find the newest matching precompiled octave file
    octfilename = '';
    for iVer = 1:numel(versionTags)
        candidate = [filename '.mexoct' versionTags{iVer} osExt bitExt];
        if exist(candidate,'file') == 2
            octfilename = candidate;
            break;
        end
    end

    fileExists = ~isempty(octfilename);

    %If it exists, we create a link with the ending "mex" which is used by octave
    if fileExists
        mexFileName = [filename '.' mexext];

        %Make the link in the right directory
        cFilePath = which(octfilename);
        [MexFolder,~,~] = fileparts(cFilePath);

        oldpath = fullfile(MexFolder,octfilename);
        linkpath = fullfile(MexFolder,mexFileName);

        %Create link

        [status,msg] = link(oldpath,linkpath);

        %We need this so the file directory cache is updated to recognize the linked file
        rehash;

        if status == 0 && exist([filename '.mex']) == 3
            matRad_cfg.dispWarning('Trying to use the precompiled mex file ''%s'' for Octave %s. This is experimental!',octfilename,ver);
            fileExists = true;
        else
            matRad_cfg.dispWarning('Could not link existing precompiled mex file %s to %s! Reason: %s',octfilename,mexFileName, msg);
            fileExists = false;
        end
    end
end

