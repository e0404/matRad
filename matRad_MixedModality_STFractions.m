matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 50000;
load 'TG119.mat'
%% 

% meta information for treatment plan (1) 
pln(1).numOfFractions  = 10;
pln(1).radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln(1).machine         = 'Generic';

% beam geometry settings
pln(1).propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln(1).propStf.gantryAngles    = [45 315]; % [?] ;
pln(1).propStf.couchAngles     = zeros(numel(pln(1).propStf.gantryAngles),1); % [?] ; 
pln(1).propStf.numOfBeams      = numel(pln(1).propStf.gantryAngles);
pln(1).propStf.isoCenter       = ones(pln(1).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(1).propDoseCalc.calcLET = 0;

pln(1).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(1).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(1).propOpt.spatioTemp      = 1;
pln(1).propOpt.STscenarios     = 3;
pln(1).propOpt.STfractions     = [2 3 5];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(1).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln(1).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(1).propDoseCalc.doseGrid.resolution = ct.resolution;
quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(1).bioParam = matRad_bioModel(pln(1).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(1).multScen = matRad_multScen(ct,scenGenType);
% 
% meta information for treatment plan (2) 
pln(2).numOfFractions  = 25;
pln(2).radiationMode   = 'photons';           % either photons / protons / helium / carbon
pln(2).machine         = 'Generic';

% beam geometry settings
pln(2).propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln(2).propStf.gantryAngles    = [0:90:359]; % [?] ;
pln(2).propStf.couchAngles     = zeros(numel(pln(2).propStf.gantryAngles),1);  % [?] ; 
pln(2).propStf.numOfBeams      = numel(pln(2).propStf.gantryAngles);
pln(2).propStf.isoCenter       = ones(pln(2).propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln(2).propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln(2).propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln(2).propOpt.spatioTemp      = 1;
pln(2).propOpt.STscenarios     = 2;
pln(2).propOpt.STfractions     = [20 5];             % can also do different spread of the fractions between scenes ( make sure sum(STfractions == numOfFractions)

% dose calculation settings
pln(2).propDoseCalc.doseGrid.resolution.x = 8; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.y = 8; % [mm]
pln(2).propDoseCalc.doseGrid.resolution.z = 8; % [mm]
% pln(2).propDoseCalc.doseGrid.resolution = ct.resolution;

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln(2).bioParam = matRad_bioModel(pln(2).radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln(2).multScen = matRad_multScen(ct,scenGenType);

% prepping cst 
% placing alpha/beta ratios in cst{:,6},
% different alpha beta ration for each obj of a structure  
sparecst = 0;

cst = matRad_prepCst(cst, sparecst);
%% Plan Wrapper
plnJO = matRad_plnWrapper(pln);
%% Stf Wrapper
stf = matRad_stfWrapper(ct,cst,plnJO);
%% Dij Calculation
 
dij = matRad_calcCombiDose(ct,stf,plnJO,cst,false);
%% Fluence optimization 
resultGUI = matRad_fluenceOptimizationJO(dij,cst,plnJO);


%% Visualization complete
%First proton vs photon vs total
%The all proton scenarios
%Then all photon scenarios

slice = 59;
qtOpt = 'physicalDose';
modalities = {'proton', 'photon'};
plan_JO = resultGUI{1}.(qtOpt) + resultGUI{2}.(qtOpt); %pln(1).numOfFractions * resultGUI{1}.(qtOpt) + pln(2).numOfFractions * resultGUI{2}.(qtOpt);

if any(plnJO.propOpt.spatioTemp)
    f = figure;
    protonScenarios = plnJO.propOpt.STscenarios(1);
    if protonScenarios>1
        for scenarioIdx = 1:protonScenarios
            subplot(3,max([plnJO.propOpt.STscenarios])+1,scenarioIdx);
            imagesc(resultGUI{1}.([qtOpt, '_STscenario_', num2str(scenarioIdx)])(:,:,slice));
            matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
            title(['proton plan scenario ', num2str(scenarioIdx)]);
            colorbar();
        end
    end
    subplot(3,max([plnJO.propOpt.STscenarios])+1,max([plnJO.propOpt.STscenarios])+1);
    imagesc(resultGUI{1}.(qtOpt)(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('total proton plan');
    colorbar();
   
    photonScenarios = plnJO.propOpt.STscenarios(2);
    if photonScenarios>1
        for scenarioIdx = 1:photonScenarios
            subplot(3,max([plnJO.propOpt.STscenarios])+1,max([plnJO.propOpt.STscenarios])+1 + scenarioIdx);
            imagesc(resultGUI{2}.([qtOpt, '_STscenario_', num2str(scenarioIdx)])(:,:,slice));
            matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
            title(['photon plan scenario ', num2str(scenarioIdx)]);
            colorbar();
        end
    end
    subplot(3,max([plnJO.propOpt.STscenarios])+1,2*(max([plnJO.propOpt.STscenarios])+1));
    imagesc(resultGUI{2}.(qtOpt)(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('total photon plan');
    colorbar();
    subplot(3,max([plnJO.propOpt.STscenarios])+1,3*(max([plnJO.propOpt.STscenarios])+1));
    imagesc(plan_JO(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
    title('total plan');
    colorbar();
end


%%% Plot photon-proton single scenarios total plan
if (protonScenarios + photonScenarios) ==2
    f = figure;

    for modalityIdx =1:2
        subplot(1,3,modalityIdx);
        imagesc(resultGUI{modalityIdx}.(qtOpt)(:,:,slice));
        matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
        title([ modalities{modalityIdx}, ' plan']);
    end

    subplot(1,3,3);
    imagesc(plan_JO(:,:,slice));
    matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

    title('Total plan');
end


%% Visulaization

slice = 59;
qtOpt = 'physicalDose';
f = figure;
plan_JO = resultGUI{1}.(qtOpt) + resultGUI{2}.(qtOpt); %pln(1).numOfFractions * resultGUI{1}.(qtOpt) + pln(2).numOfFractions * resultGUI{2}.(qtOpt);

plan_proton_total = resultGUI{1}.(qtOpt);
plan_photon        = resultGUI{2}.(qtOpt);

plan_proton_1 = resultGUI{1}.([qtOpt, '_STscenario_1']);
plan_proton_2 = resultGUI{1}.([qtOpt, '_STscenario_2']);


subplot(1,3,1);
imagesc(plan_proton_total(:,:,slice));
title('proton plan');
matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

subplot(1,3,2);
imagesc(plan_photon(:,:,slice));
title('photon plan');
matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

subplot(1,3,3);
imagesc(plan_JO(:,:,slice));
title('total plan');
matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);


f = figure;
subplot(1,4,1);
imagesc(plan_proton_total(:,:,slice));
title('total prioton plan');
matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

subplot(1,4,2);
imagesc(plan_proton_1(:,:,slice));
title('STscenario 1 plan');
matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
colorbar;
subplot(1,4,3);
imagesc(plan_proton_2(:,:,slice));
title('Stscenario 2 plan');
matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);

colorbar;
subplot(1,4,4);
imagesc(plan_proton_1(:,:,slice)./dij.STfractions{1}(1) - plan_proton_2(:,:,slice)./dij.STfractions{1}(2));
title('differences');
matRad_plotVoiContourSlice(gca(f), cst,ct.cube, 1, 1,3,slice);
colorbar;