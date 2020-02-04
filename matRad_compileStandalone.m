%Create Standalone with and without runtime

if ismac
    suffix = 'Mac';
    installSuffix = 'dmg';
    % Mac platform
elseif isunix
    suffix = 'Linux';
    installSuffix = 'install';
    % Linux platform
elseif ispc
    suffix = 'Windows';
    installSuffix = 'exe';
    % Windows platform
else
    error('Platform not supported')
end



FileName=['matRad' suffix '.prj'];

xDoc = xmlread(FileName);

%Get version and Filenames
version= replace(xDoc.getElementsByTagName('param.version').item(0).getFirstChild.getData().toCharArray()','.','_');
filename_webRT = xDoc.getElementsByTagName('param.web.mcr.name').item(0).getFirstChild.getData().toCharArray()';
filename_withRT = xDoc.getElementsByTagName('param.package.mcr.name').item(0).getFirstChild.getData().toCharArray()';

% add the version to the file name if it's different from the version
%without runtime
split_name=split(filename_webRT,'_v_');

if(size(split_name,1)==2)
    if  string(split_name(2,1)) ~= version
        filename_webRT=[string(split_name(1,1)) '_v_' version];
        xDoc.getElementsByTagName('param.web.mcr.name').item(0).getFirstChild.setData(filename_webRT);
    end
    
else
    %add version
    filename_webRT=[filename_webRT '_v_' version];
    
    xDoc.getElementsByTagName('param.web.mcr.name').item(0).getFirstChild.setData(filename_webRT);

end

%with runtime
split_name=split(filename_withRT,'_v_');

if(size(split_name,1)==2)
    if  string(split_name(2,1)) ~= version
        filename_withRT=[string(split_name(1,1)) '_v_' version];
        xDoc.getElementsByTagName('param.package.mcr.name').item(0).getFirstChild.setData(filename_withRT);
    end
    
else
    %add version
    filename_withRT=[filename_withRT '_v_' version];
    xDoc.getElementsByTagName('param.package.mcr.name').item(0).getFirstChild.setData(filename_withRT);

end



%check if the files exist, and delete them
if exist(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix], 'file') == 2
  delete(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix]);
end

if exist(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix], 'file') == 2
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
  if exist(['standalone' filesep 'for_redistribution' filesep  filename_withRT '.' installSuffix], 'file') == 2
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
  if exist(['standalone' filesep 'for_redistribution' filesep  filename_webRT '.' installSuffix], 'file') == 2
      worked = true;
      break
  end
end

if worked == false
    error('Packaging Failed!');
end
disp('Packaged compiler without runtime');


disp('Packaged both compilers');


