function fileExists = matRad_checkMexFileExists(filename) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  
fileExists = (exist(filename,'file') == 3);

if ~fileExists && strcmp(matRad_getEnvironment(),'OCTAVE')
    
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

    %Check if file exists
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

