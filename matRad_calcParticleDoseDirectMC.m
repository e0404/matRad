function resultGUI = matRad_calcParticleDoseDirectMC(MCengine,ct,cst,stf,pln,w,totalNumOfParticles)

if ~exist('totalNumOfParticles')
    totalNumOfParticles = 1e7;
end

totalPlanParticles = sum(w)*1e6;

if strcmp(MCengine, 'TOPAS')
    TopasConfig.environment = 'wsl';
    TopasConfig.command = ['wsl export TOPAS_G4_DATA_DIR=~/G4Data; ~/topas/bin/topas'];
    TopasConfig.outputType = 'binary';
    TopasConfig.fracHistories = totalNumOfParticles / totalPlanParticles;
    TopasConfig.beamProfile = 'biGaussian';
    TopasConfig.useOrigBaseData = false;
    [resultGUI, TopasConfig] = matRad_calcParticleDoseMCtopas(ct,stf,pln,w,TopasConfig);
elseif strcmp(MCengine, 'MCsquare')
    resultGUI = matRad_calcDoseDirectMCsquare(ct,stf,pln,cst,w,totalNumOfParticles);
else
    error('Error: Select a valid MonteCarlo engine!');
end
