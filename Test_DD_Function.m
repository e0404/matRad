matRad_rc
load("NEW-BOXPHANTOM-Overlap.mat")

cst{1,6}{2} = struct(mLETDoseObjectives.matRad_SquaredOverdosingmLETDose(100,120));
cst{2,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(800,60));

pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';
pln.numOfFractions  = 30;


% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
% pln(1).propStf.gantryAngles    = [ 45 0 -45]; % [?] ;
pln.propStf.gantryAngles    = [90];
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1); % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propDoseCalc.calcLET = 1;

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.spatioTemp      = 0;
% plnN(1).propOpt.STscenarios     = 2;
%pln(1).propOpt.STfractions     = [ 4 4 6 8 8];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'effect';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);
% 
% sparecst = 0;
% 
% cst = matRad_prepCst(cst, sparecst);

stf = matRad_generateStf(ct,cst,pln);
% Dij Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);
% Dirty Dose Calculation
dij = matRad_calcDirtyDose(2,dij);
result = matRad_fluenceOptimization(dij,cst,pln);
