%% Test neutron dose calculation using MCNP dose engine - modification of example 1
clear
%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
load('Patient_FRM2_submGlTumor.mat');
% now we have ct data and cst data for a new phantom
display(ct);
display(cst);

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.



%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use photons for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'photons_Generic.mat'; consequently the machine has to be set to 'Generic'
pln.radiationMode = 'neutrons';
pln.machine       = 'MCNPneutronField'; 

% MCNP parameter
pln.propMCNP.normalizationFactor = 1e5;
% Define properties for MCNP simulationr
pln.propMCNP.wantWholeZ = true;
pln.propMCNP.wantToResize = false;
pln.propMCNP.densityAir = 0.123e-3; % density air in g/cm^3 
pln.propMCNP.numberParticles = 1e3; %5e9;
% pln.propMCNP.phyMode = 'N P E H'; %'N', 'N P E H D A #';

pln.propMCNP.tallySpecifier = 'TotalDose_TMESH'; % 'KERMA_F4'; 
pln.propMCNP.tallyKeyword = 'TOTAL'; %Only needed for TMESH
% pln.propOpt.bioOptimization = 'RBExSecPartDose_MCDS_RMFmodel';
pln.propOpt.defaultLQmodel.ratioAlphaBeta = 5;

%%
% Define the biological optimization model for treatment planning along
% with the quantity that should be used for optimization. Possible model values 
% are:
% 'none':     physical optimization;
% 'constRBE': constant RBE of 1.1; 
% 'MCN':      McNamara-variable RBE model for protons; 
% 'WED':      Wedenberg-variable RBE model for protons
% 'LEM':      Local Effect Model 
% and possible quantityOpt are 'physicalDose', 'effect' or 'RBExD'.
modelName    = 'none'; %'none' 'constRBE'
quantityOpt  = 'physicalDose'; %'RBExD'; %'physicalDose';                                             

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 1;
pln.propStf.gantryAngles  =270;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 200; %5; %90;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
%pln.bioParam.model = 'constRBE';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); 

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]


%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);
%cst{2,6}{1,1}.parameters{1,1} = 1.5;

%% Dose Calculation
[dij,ct,stf,pln,cst] = matRad_calcNeutronDoseMCNP(ct,stf,pln,cst);

%% Export dij matrix
%matRad_exportDij('dij.bin',dij,stf);

%% Inverse Optimization for intensity-modulated photon therapy
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI

% %% Plot the resulting dose slice
% plane      = 3;
% slice      = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
% doseWindow = [0 max([resultGUI.physicalDose(:)])];
% 
% figure,title('phantom plan')
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);
% 
% 
% %% 
% % We export the the created phantom & dose as dicom. This is handled by the 
% % class matRad_DicomExporter. When no arguments are given, the exporter searches
% % the workspace itself for matRad-structures. The output directory can be set by
% % the property dicomDir. While the different DICOM datasets (ct, RTStruct, etc) 
% % can be exported individually, we call the wrapper to do all possible exports.
% dcmExport = matRad_DicomExporter();
% dcmExport.dicomDir = [pwd filesep 'dicomExport'];
% dcmExport.matRad_exportDicom();

