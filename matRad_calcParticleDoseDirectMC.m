function resultGUI = matRad_calcParticleDoseDirectMC(ct,stf,pln,cst,w,nCasePerBixel,MCsettings)

if exist('MCsettings')
    MCsettings = matRad_MCinit(pln,MCsettings);
else
    MCsettings = matRad_MCinit(pln,[]);
end

if strcmp(MCsettings.MCengine, 'TOPAS')
    resultGUI = matRad_calcParticleDoseMCtopas(ct,stf,pln,w,MCsettings);
elseif strcmp(MCsettings.MCengine, 'MCsquare')
    resultGUI = matRad_calcDoseDirectMCsquare(ct,stf,pln,cst,w,nCasePerBixel);
else
    error('Error: Select a valid MonteCarlo engine!');
end
