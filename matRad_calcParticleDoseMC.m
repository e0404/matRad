function dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,w,nCasePerBixel,calcDoseDirect,MCsettings)

if exist('MCsettings')
    MCsettings = matRad_MCinit(pln,MCsettings);
else
    MCsettings = matRad_MCinit(pln,[]);
end

if strcmp(MCsettings.MCengine, 'TOPAS')
    dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,w,MCsettings);
elseif strcmp(MCsettings.MCengine, 'MCsquare')
    dij = matRad_calcParticleDoseMCsquare(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
else
    error('Error: Select a valid MonteCarlo engine!');
end
