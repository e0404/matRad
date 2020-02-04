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

%Get Filenames
filename_webRT = xDoc.getElementsByTagName('param.web.mcr.name').item(0).getFirstChild.getData().toCharArray()';
filename_withRT = xDoc.getElementsByTagName('param.package.mcr.name').item(0).getFirstChild.getData().toCharArray()';

%compile with runtime
xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('false');
xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('true');

xmlwrite(FileName,xDoc);

clc
applicationCompiler('-package',FileName);

cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
% loop until 'Package finished' or 'Package failed' is found in the command window
tic
worked = false;
while toc<500
  pause ( 2 )
  %myTxt = cmdWinDoc.getText(cmdWinDoc.getStartPosition.getOffset,cmdWinDoc.getLength);
  %
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

clc
applicationCompiler('-package',FileName);

%cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
% loop until 'Package finished' or 'Package failed' is found in the command window
tic
worked = false;
while toc<500
  pause ( 2 )
  %myTxt = cmdWinDoc.getText(cmdWinDoc.getStartPosition.getOffset,cmdWinDoc.getLength);
  %
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


