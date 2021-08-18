function matRad_stopAllOpenstack()
[usr,pw,ip] = getOpenStackCredentials();

call = strjoin(['wsl plink -pw ''',pw,''' ',usr,'@',ip,' ''squeue'''],'');
[~,result] = system(call);
jobList = strsplit(result,'\n')';
jobIDs = cell2mat(cellfun(@(v)str2double(v),regexp(jobList,'\d*','Match'),'UniformOutput',false));
jobIDs = jobIDs(:,1);

call = strjoin(['wsl plink -pw ''',pw,''' ',usr,'@',ip,' ''scancel ',num2str(jobIDs'),''''],'');
system(call);

end

