function [versionString, matRadVer] = matRad_version(matRadRoot)
% matRad function to get the current matRad version
% (and git information when used from within a repository
%
% call
%   [versionString,matRadVer] = matRad_version()
%   [versionString,matRadVer] = matRad_version(matRadRoot)
%
% input
%   matRadRoot:     Optional Root Directory. This is for call in matRad
%                   initialization when MatRad_Config is not yet available
% output
%   versionString:  Readable string build from version information
%   matRadVer:      struct with version information
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hardcoded version name / numbers
matRadVer.name = 'Cleve';
matRadVer.major = 3;
matRadVer.minor = 2;
matRadVer.patch = 0;
matRadVer.revision = 65534; % Fallback value for unknown revision (e.g., when not running from a git repository)
matRadVer.branch = [];
matRadVer.commitID = [];

% Determine repo root
if nargin == 1
    repoDir = matRadRoot;
else
    try
        matRad_cfg = MatRad_Config.instance();
        repoDir = matRad_cfg.matRadRoot;
    catch
        repoDir = fileparts(mfilename('fullpath'));
    end
end

% Retrieve branch / commit / revision via system git (preferred)
gitAvailable = false;
gitRepo = exist(fullfile(repoDir, '.git'), 'dir') == 7;

if gitRepo
    try
        git = @(cmd) system(sprintf('git -C "%s" %s', repoDir, cmd));

        [st1, branch]   = git('rev-parse --abbrev-ref HEAD');
        [st2, commitID] = git('rev-parse HEAD');
        branch   = strtrim(branch);
        commitID = strtrim(commitID);

        if st1 == 0 && st2 == 0 && numel(commitID) == 40
            gitAvailable = true;
            matRadVer.commitID = commitID;
            if strcmp(branch, 'HEAD')
                matRadVer.branch = 'DETACHED';
            else
                matRadVer.branch = branch;
            end

            % Count commits since version tag -> revision number
            tag = sprintf('v%d.%d.%d', matRadVer.major, matRadVer.minor, matRadVer.patch);
            [stRev, revCount] = git(sprintf('rev-list --count %s..HEAD', tag));
            if stRev == 0
                matRadVer.revision = str2double(strtrim(revCount));
            end
        end
    catch
        % system git not available, will attempt file-based fallback below
    end

    % Fallback: read branch / commit from .git files directly (no revision)
    if ~gitAvailable
        try
            repoGitDir = fullfile(repoDir, '.git');
            headText = fileread(fullfile(repoGitDir, 'HEAD'));

            % Test if detached head (HEAD contains 40 hex SHA1 commit ID)
            i = regexp(headText, '[a-f0-9]{40}', 'once');
            if ~isempty(i)
                matRadVer.branch   = 'DETACHED';
                matRadVer.commitID = headText(1:40);
            else % HEAD contains reference to branch
                headParse = textscan(headText, '%s');
                refHead   = headParse{1}{2};
                refParse  = strsplit(refHead, '/');
                matRadVer.branch = strjoin(refParse(3:end), '/');

                % Read commit ID from ref path
                refID = fileread(fullfile(repoGitDir, strjoin(refParse, filesep)));
                matRadVer.commitID = refID(1:40);
            end
        catch
            % Git repo information could not be read - keep defaults
        end
    end

    % revision == 0 means we are sitting exactly on the version tag (release)
    tagged = matRadVer.revision == 0;

    % Create a readable string
    if ~isempty(matRadVer.branch) && ~isempty(matRadVer.commitID) && ~tagged
        gitString = sprintf('(%s-%s)', matRadVer.branch, matRadVer.commitID(1:8));
    end
    versionString = sprintf('"%s" v%d.%d.%d.%d%s', matRadVer.name, matRadVer.major, matRadVer.minor, matRadVer.patch, matRadVer.revision, gitString);
else
    versionString = sprintf('"%s" v%d.%d.%d', matRadVer.name, matRadVer.major, matRadVer.minor, matRadVer.patch);
end
% This checks if no explicit assignment is done in which case the version is printed.
if nargout == 0
    disp(['You are running matRad ' versionString]);
end

end
