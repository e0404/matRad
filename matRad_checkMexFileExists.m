function fileExists = matRad_checkMexFileExists(filename,noLinkOctave) 
% Checks if a matching mex file exists, and creates a link if a matching 
% octave mex file is found
%
% call
%   matRad_checkMexFileExists(filename)
%
%
% References
%
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

  
%Explicitly check for matching mex file (id 3)
fileExists = (exist(filename,'file') == 3);

if nargin < 2
    noLinkOctave = false;
end

if ~fileExists && strcmp(matRad_getEnvironment(),'OCTAVE') && ~noLinkOctave
    
    [~,maxArraySize] = computer();
    
    if maxArraySize > 2^31
        bitExt = '64';
    else
        bitExt = '32';
    end
    
    if ispc
        systemext = 'mexoctw';
    elseif ismac 
        systemext = 'mexoctmac';
    elseif isunix
        systemext = 'mexocta';
    else
        error('No mex file for octave for your operating system. Compile it yourself.');
    end
    
    octfilename = [filename '.' systemext bitExt];

    %Check if matching octave file exists
    fileExists = (exist(octfilename,'file') == 2);
    
    if fileExists
        %Make the link in the right directory
        cFilePath = which(octfilename);
        [MexFolder,~,~] = fileparts(cFilePath);
        
        currFolder = pwd;
        
        cd(MexFolder);
        
        %Create link
        link(octfilename, [filename '.mex']);
        
        cd(currFolder);
        warning('You use a precompiled mex with Octave. This is experimental!');
    end

end

