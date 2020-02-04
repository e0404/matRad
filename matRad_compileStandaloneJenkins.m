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


disp('Packaged compiler without runtime');


%compile with runtime
xDoc.getElementsByTagName('param.web.mcr').item(0).getFirstChild.setData('false');
xDoc.getElementsByTagName('param.package.mcr').item(0).getFirstChild.setData('true');

xmlwrite(FileName,xDoc);

disp('Packaged both compilers');



