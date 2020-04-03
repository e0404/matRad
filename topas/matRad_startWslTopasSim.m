function matRad_startWslTopasSim(ct,stf,pln,w,wslDistribution,userName,TopasConfig)

if ~exist('wslDistribution','var')
    wslDistribution = 'Ubuntu-18.04';
end

if ~exist('userName','var')
    userName = 'homolka';
end

TopasConfig.filepath = ['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\'];

load([pln.radiationMode,'_',pln.machine]);
topasBaseData = MatRad_TopasBaseData(machine);
% generate Topas data at the topas installation directory
topasBaseData.writeTopasData(ct,stf,pln,w,TopasConfig);

startSim = ['wsl -d ',wslDistribution,' cd ~/; source setup_env.sh; cd ~/matfiles/MCexport/; startTopasSim &'];
system(startSim)

end


