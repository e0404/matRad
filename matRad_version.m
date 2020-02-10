function [versionString,matRadVer] = matRad_version()
% matRad function to get the current matRad version (and git information when 
% used from within a repository
% 
% call
%   [versionString,matRadVer] = matRad_version()
%
% input
%   
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Hardcoded version name / numbers
matRadVer.name = 'Blaise';
matRadVer.major = 2;
matRadVer.minor = 5;
matRadVer.patch = 0;

%Retreive branch / commit information from current git repo if applicable
try
    %read HEAD file to point to current ref / commit
    repoGitDir = [fileparts(mfilename('fullpath')) filesep '.git'];   
    headText = fileread([repoGitDir filesep 'HEAD']);
    
    %Test if detached head (HEAD contains 40 hex SHA1 commit ID)
    i = regexp(headText,'[a-f0-9]{40}');
    if ~isempty(i)
        matRadVer.branch = 'DETACHED';
        matRadVer.commitID = headText(1:40);        
    else %HEAD contains reference to branch
        headParse = textscan(headText,'%s');    
        refHead = headParse{1}{2};
        refParse = strsplit(refHead,'/');
        matRadVer.branch = refParse{end};

        %Read ID from ref path
        refID = fileread([repoGitDir filesep strjoin(refParse,filesep)]);
        matRadVer.commitID = refID(1:40);
    end
catch
    %Git repo information could not be read, set to empty
    matRadVer.branch = [];
    matRadVer.commitID = [];    
end

%Create a readable string
%Git path first
gitString = '';
if ~isempty(matRadVer.branch) && ~isempty(matRadVer.commitID)
    gitString = sprintf(' (%s-%s)',matRadVer.branch,matRadVer.commitID(1:8));  
end

%Full string
versionString = sprintf('v%d.%d.%d "%s"%s',matRadVer.major,matRadVer.minor,matRadVer.patch,matRadVer.name,gitString);

%This checks if no explicit assigment is done in which case the version is printed.
if nargout == 0
    disp(['You are running matRad ' versionString]);
end

end
