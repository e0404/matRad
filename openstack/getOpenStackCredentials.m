function [usr,pw,ip] = getOpenStackCredentials()

fID = fopen('access');
access = convertCharsToStrings(fread(fID,'*char'));
fclose(fID);

access = strsplit(access);
usr = access(1);
pw = access(2);
ip = access(3);

end

