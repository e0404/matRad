
matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
ctDim = [200,200,100]; % x,y,z dimensions
ctResolution = [2,2,3]; % x,y,z the same here!

builder = matRad_PhantomBuilder(ctDim,ctResolution,1);


objective3 = struct(DoseObjectives.matRad_SquaredDeviation(800,45));
objective2 = struct(DoseObjectives.matRad_SquaredOverdosing(10,0));
objective1 = struct(DoseObjectives.matRad_SquaredOverdosing(400,0));

builder.addBoxTarget('Volume3',[20,20,20],'objectives',objective3)
builder.addBoxOAR('Volume1',[60,30,40],'offset',[0 -15 0],'objectives',objective1);
builder.addBoxOAR('Volume2',[60,30,40],'offset',[0 15 0],'objectives',objective2);
builder.addBoxOAR('Volume4',[100,120,80],'offset',[0 0 0]);
%%

[ct,cst] = builder.getctcst();

cst{1,6}{2} = DoseConstraints.matRad_MinMaxDose(40,52);
%{
cst{1,6}{2} = DoseConstraints.matRad_MinMaxDose(0,40);

cst{2,6}{2} = DoseConstraints.matRad_MinMaxDose(0,40);
%}
matRadGUI
%%
%pln.radiationMode = 'photons';            
%changed
pln.radiationMode = 'protons';
pln.machine       = 'Generic';

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
modelName    = 'MCN';
quantityOpt  = 'effect';                                             

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [0,180];
%pln.propStf.gantryAngles  = [0:70:355];
pln.propStf.couchAngles   = zeros(size(pln.propStf.gantryAngles));
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
%pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*[200 200 75]
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); 

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 4; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 4; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 4; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);
 
%% Dose Calculation
dij = matRad_calcParticleDose(ct,stf,pln,cst);
%dij = matRad_calcPhotonDose(ct,stf,pln,cst);
%%
%%
%matRad_UIInterpolation(data,dij,pln,ct,cst)
%'a'
%%


%[resultGUI,optimizer] = matRad_fluenceOptimization(dij,cst,pln);
%matRadGUI;
%load('dataFull.mat')
%retStruct.modcst = cst2;

%matRad_UIInterpolation(retStruct,dij,pln,ct,cst,retStruct.optiProb)
%%
retStruct = matRad_ParetoOptimization(dij,cst,pln,20);
%%
%load('correct100.mat')
%load('data50protons1constr.mat')

%%
matRad_UIInterpolation(retStruct,dij,pln,ct,cst,retStruct.optiProb)

%%
matRad_plotParetoSurface(retStruct)
%%
ps = retStruct.finds
%%
fValsMod = matRad_generateDummyPoints(ps,ones(3,1))
%%
[k] = convhulln(fValsMod)
figure
trisurf(k,fValsMod(:,1),fValsMod(:,2),fValsMod(:,3),'FaceColor','cyan')
%%

[k] = convhulln(ps);
figure
trisurf(k,ps(:,1),ps(:,2),ps(:,3),'FaceColor','cyan')
%scatter3(ps(:,1),ps(:,2),ps(:,3),'MarkerEdgeColor','black',...
%        'MarkerFaceColor',[0 0 0])
        %scatter3(psRed(:,1),psRed(:,2),psRed(:,3),'MarkerEdgeColor','black',...
        %        'MarkerFaceColor',[1 1 0])
        %scatter3(fNew(:,1),fNew(:,2),fNew(:,3),'filled','MarkerFaceColor','red')
%%
matRadGUI