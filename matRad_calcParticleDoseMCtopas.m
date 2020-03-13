function dij = matRad_calcParticleDoseMCtopas(ct,stf,pln,w,MCsettings)




basematerial = '';
if ~exist('machine') || ~isfield(machine.meta,'basematerial')
  warning('Base material not defined in base data. Using Water')
  basematerial = 'Water';
else
  basematerial = machine.meta.basematerial;
end

matRad_exportCtTOPAS(ct, MCsettings.runsPath, basematerial)
