function matRad_writeToOpenstack(folderToCopy)
[usr,pw,ip] = getOpenStackCredentials();

[~,wslFolder] = system(['wsl wslpath -a ''',folderToCopy,'''']);
wslFolder = convertCharsToStrings(wslFolder(1:end-1));
call = strjoin(['wsl pscp -pw ''',pw,''' -r ',wslFolder,' ',usr,'@',ip,':/nfs/home/',usr,'/files/'],'');
[status,errMsg] = system(call);

if status == 1
   warning(convertCharsToStrings(errMsg(1:end-2)));
else
   disp('Successfully written to OpenStack.');
end
end

