load PROSTATE.mat

pln.radiationMode                       = 'protons';  
pln.machine                             = 'Generic';
pln.propOpt.bioOptimization             = 'none';    
pln.numOfFractions                      = 1;
pln.propStf.gantryAngles                = [90 270];
%pln.propStf.gantryAngles               = [0];
pln.propStf.couchAngles                 = zeros(1, numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth                  = 5;
pln.propStf.numOfBeams                  = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter                   = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propDoseCalc.doseGrid.resolution    = ct.resolution;

pln.propOpt.runDAO = 0;

stf = matRad_generateStf(ct,cst,pln);
dij = matRad_calcParticleDose(ct,stf,pln,cst);

load prostate_super_cst.mat

save PROSTATE_super.mat