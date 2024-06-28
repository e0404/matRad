function matRad_compileStandalone(varargin)
% Compiles the standalone exectuable & packages installer using Matlab's
% Compiler Toolbox
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


p = inputParser;

p.addParameter('isRelease',false,@islogical);               %By default we compile a snapshot of the current branch
p.addParameter('compileWithRT',false,@islogical);           %By default we don't package installers with runtime
p.addParameter('compileWithoutRT',true,@islogical);         %By default we do package installers without runtime
p.addParameter('projectFile','matRad.prj',@(x) exist(x,2)); %Default template prj file
p.addParameter('override',false,@islogical);                %If set we do not care for any argumetn but use the bare *.prj file specified

p.parse(varargin{:});

isRelease = p.Results.isRelease;
compileWithRT = p.Results.compileWithRT;
compileWithoutRT = p.Results.compileWithoutRT;
projectFile = p.Results.projectFile;
override = p.Results.override;


if ~compileWithRT && ~compileWithoutRT
    error('No target specified!');
end

[~,versionFull] = matRad_version();

%Create Standalone with and without runtime
xDoc = xmlread(projectFile);

if ~override
    
    if ismac % Mac platform
        suffix = 'Mac';
        installSuffix = 'app';
        xDoc.getElementsByTagName('unix').item(0).getFirstChild.setData('true');
        xDoc.getElementsByTagName('mac').item(0).getFirstChild.setData('true');
        xDoc.getElementsByTagName('windows').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('linux').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('arch').item(0).getFirstChild.setData('maci64');
        
    elseif isunix % Linux platform
        suffix = 'Linux';
        installSuffix = 'install';
        xDoc.getElementsByTagName('unix').item(0).getFirstChild.setData('true');
        xDoc.getElementsByTagName('mac').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('windows').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('linux').item(0).getFirstChild.setData('true');
        xDoc.getElementsByTagName('arch').item(0).getFirstChild.setData('glnxa64');
        
    elseif ispc % Windows platform
        suffix = 'Win';
        installSuffix = 'exe';
        xDoc.getElementsByTagName('unix').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('mac').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('windows').item(0).getFirstChild.setData('true');
        xDoc.getElementsByTagName('linux').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('arch').item(0).getFirstChild.setData('Win64');
    else
        error('Platform not supported')
    end
    
    %Set MATLAB root
    xDoc.getElementsByTagName('root').item(0).getFirstChild.setData(matlabroot);
    
    %Replace the file separator in the file tags
    files = xDoc.getElementsByTagName('file');
    for i=0:files.getLength - 1
        fileName=strrep(strrep(files.item(i).getFirstChild.getData().toCharArray()','/',filesep),'\',filesep);
        files.item(i).getFirstChild.setData(fileName);
        
        %change the file separator in the 'location' attribute of the file tag
        %(if it exists)
        if hasAttributes(files.item(i))
            location=strrep(strrep(files.item(i).getAttribute('location').toCharArray()','/',filesep),'\',filesep);
            files.item(i).setAttribute('location',location);
        end
    end
    
    %Replace the file separator in the attributes and children tags of configuration
    config = xDoc.getElementsByTagName('configuration');
    config.item(0).setAttribute('file',['${PROJECT_ROOT}' filesep projectFile]);
    config.item(0).setAttribute('location','${PROJECT_ROOT}');
    config.item(0).setAttribute('preferred-package-location',['${PROJECT_ROOT}' filesep 'standalone' filesep 'for_redistribution']);
    
    tags = config.item(0).getChildNodes;
    
    for i=0:tags.getLength - 1
        if ~isempty(tags.item(i).getFirstChild)
            fileName=strrep(strrep(tags.item(i).getFirstChild.getData().toCharArray()','/',filesep),'\',filesep);
            tags.item(i).getFirstChild.setData(java.lang.String(fileName));
        end
    end
    
    %Get version and Filenames
    %version= replace(xDoc.getElementsByTagName('param.version').item(0).getFirstChild.getData().toCharArray()','.','_');
    
    vernum = sprintf('%d.%d.%d',versionFull.major,versionFull.minor,versionFull.patch);
    
    tmp_verstr = ['v' vernum];
    if ~isRelease && ~isempty(versionFull.commitID) && ~isempty(versionFull.branch)
        branchinfo = [versionFull.branch '-' versionFull.commitID(1:8)];
        tmp_verstr = [tmp_verstr '_' branchinfo];
        addDescription = ['!!!matRad development version build from branch/commit ' branchinfo '!!!'];
        desc = xDoc.getElementsByTagName('param.description').item(0).getFirstChild.getData().toCharArray()';
        desc = [addDescription ' ' desc];
        xDoc.getElementsByTagName('param.description').item(0).getFirstChild.setData(desc);
    end
    
    %The Application compiler only allows major and minor version number
    vernum_appCompiler = sprintf('%d.%d',versionFull.major,versionFull.minor);
    
    filename_webRT =['matRad_installer' suffix '64_' tmp_verstr];
    filename_withRT =['matRad_installer' suffix '64_wRT_' tmp_verstr];
    xDoc.getElementsByTagName('param.web.mcr.name').item(0).getFirstChild.setData(filename_webRT);
    xDoc.getElementsByTagName('param.package.mcr.name').item(0).getFirstChild.setData(filename_withRT);
    xDoc.getElementsByTagName('param.version').item(0).getFirstChild.setData(vernum_appCompiler);
    
    
    
    %check if the files exist, and delete them
    if isfolder(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix])
        rmdir(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix],'s')
        rehash
    end
    
    if isfile(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix])
        delete(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix]);
    end
    
    if isfolder(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix])
        rmdir(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix],'s')
        rehash
    end
    if isfile(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix])
        delete(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix]);
    end
    
    %compile with runtime
    if compileWithRT
        xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('false');
        xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('true');
        
        xmlwrite(projectFile,xDoc);
        
        applicationCompiler('-package',projectFile);
        
        % loop until the file is created
        waitForCompilation(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix],1000);
    end
    
    if compileWithoutRT
        %compile without runtime
        xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('true');
        xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('false');
        
        xmlwrite(projectFile,xDoc);
        
        applicationCompiler('-package',projectFile);
        
        % loop until the file is created
        waitForCompilation(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix],1000);
    end
    disp('Done');
else %No automatic compilation from template, override with prj file
    applicationCompiler('-package',projectFile);   
    %No waiting here because we do not know the output file
end

end

function waitForCompilation(path,waittime)
tic
worked = false;
while toc<waittime
    pause(2);
    if isfolder(path) || isfile(path)
        worked = true;
        break
    end
end

if worked == false
    error('Packaging Failed!');
end
disp('Packaged compiler without runtime');
end