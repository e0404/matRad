function TopasConfig = matRad_startWslTopasSim(TopasConfig)

if ~isfield('TopasConfig','wslDistribution')
    TopasConfig.wslDistribution = 'Ubuntu-18.04';
end

if ~isfield('TopasConfig','wslUserName')
    TopasConfig.wslUserName = 'homolka';
end

if ~isfield('TopasConfig','wslTopasFolder')
    TopasConfig.wslTopasFolder = '/matfiles/MCexport';
end

TopasConfig.topasFilepath = ['\\wsl$\',TopasConfig.wslDistribution,'\home\',TopasConfig.wslUserName,TopasConfig.wslTopasFolder];

copyfile(TopasConfig.matRadFolder,TopasConfig.topasFilepath)


%startSim = ['wsl -d ',wslDistribution,' cd ~/; source setup_env.sh; cd ~',topasFolder,'/; startTopasSim'];
startSim = ['wsl -d ',TopasConfig.wslDistribution,' cd ~/; source setup_env.sh; cd ~',TopasConfig.wslTopasFolder,'/; ../../topas/topas matRad_plan_field1_run1.txt'];
system(startSim)

end


