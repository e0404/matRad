%Create Standalone with and without runtime

FileName='matRad.prj';
xDoc = xmlread(FileName);

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

%FileName=['matRad' suffix '.prj'];

%Get version and Filenames
version= replace(xDoc.getElementsByTagName('param.version').item(0).getFirstChild.getData().toCharArray()','.','_');
filename_webRT =['matRad_installer' suffix '64_v_' version];
filename_withRT =['matRad_installer' suffix '64_wRT_v_' version];
xDoc.getElementsByTagName('param.web.mcr.name').item(0).getFirstChild.setData(filename_webRT);
xDoc.getElementsByTagName('param.package.mcr.name').item(0).getFirstChild.setData(filename_withRT);
%filename_webRT = xDoc.getElementsByTagName('param.web.mcr.name').item(0).getFirstChild.getData().toCharArray()';
%filename_withRT = xDoc.getElementsByTagName('param.package.mcr.name').item(0).getFirstChild.getData().toCharArray()';


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
xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('false');
xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('true');

xmlwrite(FileName,xDoc);

applicationCompiler('-package',FileName);

% loop until the file is created
tic
worked = false;
while toc<500
  pause ( 2 )
  if isfolder(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix]) || ...
      isfile(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix])
      worked = true;
      break
  end
end
if worked == false
    error('Packaging Failed!');
end
disp('Packaged compiler with runtime');

%compile without runtime
xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('true');
xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('false');

xmlwrite(FileName,xDoc);

applicationCompiler('-package',FileName);

% loop until the file is created
tic
worked = false;
while toc<500
  pause ( 2 )
  if isfolder(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix]) || ...
      isfile(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix])
      worked = true;
      break
  end
end

if worked == false
    error('Packaging Failed!');
end
disp('Packaged compiler without runtime');


disp('Packaged both compilers');


