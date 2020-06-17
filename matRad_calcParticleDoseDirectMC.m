function resultGUI = matRad_calcParticleDoseDirectMC(MCengine,ct,cst,stf,pln,w,totalNumOfParticles)

if ~exist('totalNumOfParticles')
    totalNumOfParticles = 1e7;
end

if strcmp(MCengine, 'TOPAS')
    TopasConfig.environment = 'wsl';
    [resultGUI, TopasConfig] = matRad_calcParticleDoseMCtopas(ct,stf,pln,w,TopasConfig);
elseif strcmp(MCengine, 'MCsquare')
    resultGUI = matRad_calcDoseDirectMCsquare(ct,stf,pln,cst,w,totalNumOfParticles);
else
    error('Error: Select a valid MonteCarlo engine!');
end
