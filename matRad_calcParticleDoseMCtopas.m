function [topasDose] = matRad_calcParticleDoseMCtopas(ct,stf,pln,w,wslDistribution,userName,TopasConfig)

%% sending data to topas

waitForSystem = 'false';

if ~exist('wslDistribution','var')
    wslDistribution = 'Ubuntu-18.04';
end

if ~exist('userName','var')
    userName = 'homolka';
end

TopasConfig.filepath = ['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\'];

load([pln.radiationMode,'_',pln.machine]);
topasBaseData = MatRad_TopasBaseData(machine,TopasConfig);
% generate Topas data at the topas installation directory
topasBaseData.writeTopasData(ct,stf,pln,w);

files = dir(fullfile(TopasConfig.filepath,'*.txt'));

paramFiles = [];
for k = 3:length(files)
    paramFiles = [paramFiles,'; ../../topas/topas ',files(k).name];
end
if waitForSystem
startSim = ['wsl -d ',wslDistribution,' cd ~/; source setup_env.sh; cd ~/matfiles/MCexport/',paramFiles];    

% receiving data from topas
topasDose = matRad_getWslTopasData(wslDistribution,userName,TopasConfig);
else
startSim = ['wsl -d ',wslDistribution,' cd ~/; source setup_env.sh; cd ~/matfiles/MCexport/',paramFiles,' &'];    
end

system(startSim)



end


