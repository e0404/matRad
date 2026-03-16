function matRad_cleanRepFromRUNTPE(inputExtension) 
% matRad_cleanRepFromRUNTPE: function cleans up directory by deleting MCNP
% files with extension given by input.
% Delete RUNTPE by default.

if nargin < 1
    inputExtension = 'r';
    filelist = dir(strcat('MCNPrunfile_*bixel', inputExtension));
else
    filelist = dir(strcat('MCNPrunfile_*bixel', inputExtension));
end

if isempty(filelist)
    warning(strcat('No files: ', ' MCNPrunfile_*bixel', inputExtension, ' were found in current directory.'))
else
    delete(filelist.name)
end
end