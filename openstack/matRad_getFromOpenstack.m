function matRad_getFromOpenstack(folderToGet,folderToWrite)
[usr,pw,ip] = getOpenStackCredentials();

[~,wslFolder] = system(['wsl wslpath -a ''',folderToWrite,'''']);
wslFolder = wslFolder(1:end-1);

% make folder if not already existent
if ~isfolder(folderToWrite)
    system(['wsl mkdir -p ',wslFolder]);
end

call = strjoin(['wsl pscp -pw ''',pw,''' -r ',usr,'@',ip,':/nfs/home/',usr,'/files/',folderToGet,' ',wslFolder],'');
[status,errMsg] = system(call);

if status == 0
    disp('Successfully copied from OpenStack.');
else
    warning(convertCharsToStrings(errMsg(1:end-2)));
end
end

