%Create Standalone with and without runtime

if ismac
    suffix = 'Mac';
    % Mac platform
elseif isunix
    suffix = 'Linux';
    % Linux platform
elseif ispc
    suffix = 'Windows';
    % Windows platform
else
    error('Platform not supported')
end



FileName=['matRad' suffix '.prj'];
xDoc = xmlread(FileName);

%compile without runtime
xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('true');
xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('false');

xmlwrite(FileName,xDoc);

clc
applicationCompiler('-package',FileName);

cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
% loop until 'Package finished' or 'Package failed' is found in the command window
tic
while toc<500
  pause ( 2 )
  myTxt = cmdWinDoc.getText(cmdWinDoc.getStartPosition.getOffset,cmdWinDoc.getLength);
  %
  if ~isempty ( strfind ( myTxt, 'Package finished' ) )
    break
  end
  if ~isempty ( strfind ( myTxt, 'Package failed' ) )
    break
  end
end

disp('Packaged compiler without runtime');


%compile with runtime
xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('false');
xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('true');

xmlwrite(FileName,xDoc);

clc
applicationCompiler('-package',FileName);

cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
% loop until 'Package finished' or 'Package failed' is found in the command window
tic
while toc<500
  pause ( 2 )
  myTxt = cmdWinDoc.getText(cmdWinDoc.getStartPosition.getOffset,cmdWinDoc.getLength);
  %
  if ~isempty ( strfind ( myTxt, 'Package finished' ) )
    break
  end
  if ~isempty ( strfind ( myTxt, 'Package failed' ) )
    break
  end
end

disp('Packaged both compilers');



