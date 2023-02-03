function fileExists = matRad_checkMexFileExists(filename,linkOctave)
% Checks if a matching mex file exists, and can create a link
% if a matching custom, system specific precompiled octave mex file is found
%
% call
%   matRad_checkMexFileExists(filename)
%   matRad_checkMexFileExists(filename,linkOctave)
%
% input:
%   filename:   name of the mex file (without extension)
%   linkOctave: (optional: default true) If set to true, the function will check
%               for a custom build mex file for octave with our custom extension
%               mexoct<version><system>. If such a file exists, it will create
%               a link to the file as filename.mex since octave not uses system-
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
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
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
fileExists = (exist(filename,'file') == 3);

%For octave we have experimental precompiled files for Octave 5 64 bit
[env,ver] = matRad_getEnvironment();

if ~fileExists && matRad_cfg.isOctave && linkOctave

    %versionstring
    verStripped = erase(ver,'.');

    %Check Architecture
    [~,maxArraySize] = computer();

    if maxArraySize > 2^31
        bitExt = '64';
    else
        bitExt = '32';
    end


    systemext = ['mexoct' verStripped];
    if ispc
        systemext = [systemext 'w'];
    elseif ismac
        systemext = [systemext 'mac'];
    elseif isunix
        systemext = [systemext 'a'];
    else
        %No file for unknown operating system
        fileExists = false;
        return;
    end

    %Build our octfilename
    octfilename = [filename '.' systemext bitExt];

    %Check if matching precompiled octave file exists
    fileExists = (exist(octfilename,'file') == 2);

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
            matRad_cfg.dispWarning('Trying to use a precompiled mex for Octave %s. This is experimental!',ver);
            fileExists = true;
        else
            matRad_cfg.dispWarning('Could not link existing precompiled mex file %s to %s! Reason: %s',octfilename,mexFileName, msg);
            fileExists = false;
        end
    end
end

