load TG119.mat

pln.radiationMode                       = 'photons';  
pln.machine                             = 'Generic';
  
pln.numOfFractions                      = 1;
pln.propStf.gantryAngles                = [0:72359];
pln.propStf.couchAngles                 = zeros(1, numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth                  = 5;
pln.propStf.numOfBeams                  = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter                   = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propDoseCalc.doseGrid.resolution    = ct.resolution;

pln.propOpt.runDAO = 0;

stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcPhotonDose(ct,stf,pln,cst);


save TG119_super.mat