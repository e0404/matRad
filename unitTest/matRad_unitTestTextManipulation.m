function matRad_unitTestTextManipulation(filename, string1, string2, path)
% matRad function for manipulating certain parameters in matRad pipeline
% so that it would be less computational expensive when it's being 
% used in continuous integration builds.
%
% call
%   matRad_unitTestTextManipulation(filename, string1, string2, path)
%
% input
%     filename: file(s) to be manipulated, pass multiple files as cell array
%     string1: string to be changed
%     string2: string to be replaced
%     path: (optional) path where the function is located. default is the
%     parent directory '../'
%
%
% References
%   -   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    path = ['..' filesep];
end

if ~iscell(filename)
    filename = {filename};
end

for fIx = 1:numel(filename)
    
    currFilename = filename{fIx};
    
    fid=fopen([path currFilename]);
    fo=fopen('tempFile.m','w');
    tline = fgetl(fid);
    
    while ischar(tline)
        
        if strfind(tline, string1)
            fprintf(fo, '%s\n', string2);
        else
            fprintf(fo, '%s\n', tline);
        end
        tline = fgetl(fid);
    end
    
    
    fclose(fid);
    fclose(fo);
    
    
    movefile('tempFile.m', [path currFilename], 'f');
end
end

