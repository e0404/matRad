function [ct, cst, pln, dij, resultGUI, DoseDis, D95, HI, accDose_Opt, accDose_WC] = ConvVsWC()

%% Calculation of a 3D WC Plan and conventional IMPT Plan and comparison of the effect of organ motion on both plans


load('Liver_DS221.mat')

% check that PTV is target objective in cst (not GTV or ITV)
%%
% meta information for treatment plan
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'HITfixedBL';

% beam geometry settings
pln.propStf.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longSpotSpacing = 5;      % only relevant for HIT machine, not generic
pln.propStf.gantryAngles    = [240 320]; 
pln.propStf.couchAngles     = [0 0]; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

%optimization settings
pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

quantityOpt  = 'RBExD';     % options: physicalDose, effect, RBExD
modelName    = 'constRBE';             % none: for photons, protons, carbon            % constRBE: constant RBE 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions

%scenGenType  = 'nomScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'       


% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);


for i=1:size(cst,1)
    if isfield(cst{i,6}, 'robustness')
        cst{i,6}.robustness ='VWWC';
    end
end

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'wcScen');

%[cst,pln] = matRad_setPlanUncertainties(ct,cst,pln);
     
pln.multScen.shiftSize = [5 5 5];        % maximum shift [mm]  % prostate cases 5mm otherwise 3mm
pln.multScen = pln.multScen.matRad_createValidInstance;


%% generate steering file
%steering file is generated for PTV, WC and 10 scenarios
stf = matRad_generateStf(ct,cst,pln);

%% hack to only do dose calculation for first dose ct
ct.numOfCtScen = 1;
pln.multScen.numOfCtScen = 1;

%% WC Opt nur für ITV, conv Opt auf PTV
disp('check that correct objectives are defined in cst. Attention if target names are not ITV and PTV')
for i=1:size(cst,1)
        if strcmp(cst{i,2},'PTV')
            index_PTV = i;
        end
        if strcmp(cst{i,2},'ITV')
            index_ITV = i; 
        end
end

cst{index_ITV,6} = cst{index_PTV,6};
cst{index_PTV,6}= [];

%%
dij = matRad_calcParticleDose(ct,stf,pln,cst,false);

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% post processing
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf);   

%% export Plan
matRad_export_HITXMLPlan_modified('Liver_WC240320_bf_5mm_new',  pln, stf, resultGUI, 'backforth') 

weights_WC = resultGUI.w;
DoseDis.WC240320 = resultGUI.RBExD;

%%
ct.numOfCtScen = 10;
pln.multScen.numOfCtScen = 10;

%% WC Opt nur für ITV, conv Opt auf PTV
cst{index_PTV,6} = cst{index_ITV,6};
cst{index_ITV,6}= [];

%% set plan uncertainties for robust optimization

for i=1:size(cst,1)
    if isfield(cst{i,6}, 'robustness')
        cst{i,6}.robustness ='none';
    end
end

pln.multScen = matRad_multScen(ct,'nomScen');

%%
clear dij;
dij = matRad_calcParticleDose(ct,stf,pln,cst,false);

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% post processing
resultGUI = matRad_postprocessing(resultGUI, dij, pln, cst, stf);   

weights_Opt = resultGUI.w;
DoseDis.Opt240320 = resultGUI.RBExD;

%% export Plan
matRad_export_HITXMLPlan_modified('Liver_240320_bf_5mm_new',  pln, stf, resultGUI, 'backforth') 

%% calc 4D dose
resultGUI.w = weights_Opt;
resultGUI.RBExD = DoseDis.Opt240320;

[resultGUI, delivery] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, 'Liver_240320_bf_5mm_new');

DoseDis.M240320 = resultGUI.accRBExD;

%% motion for WC plan

resultGUI.w = weights_WC;
resultGUI.RBExD = DoseDis.WC240320;

[resultGUI, delivery] = matRad_calc4dDose(ct, pln, dij, stf, cst, resultGUI, 'Liver_WC240320_bf_5mm_new');

DoseDis.MWC240320 = resultGUI.accRBExD;

%% Dose statistics
for i=1:size(cst,1)
        if strcmp(cst{i,2},'pseudoGTV')
            indices     = cst{i,4}{1}; 
        end
end

 refVol = [5 95]; %5% und 95%
 numOfVoxels = numel(indices);
   
 % get Dose, dose is sorted to simplify calculations
    relevantDose = DoseDis.Opt240320;
    doseInVoi    = sort(DoseDis.Opt240320(indices));
    
    %refGy = round([0.4 0.6 0.8] * max(relevantDose(:)) * 10)/10;
    refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
      
    for runDX = 1:numel(refVol)
       QI.(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));          
    end    
    
    D95.Opt240320 = QI.D95;        
    HI.Opt240320 = (DX(2)/DX(98));
    
     % get Dose, dose is sorted to simplify calculations
    relevantDose = DoseDis.M240320;
    doseInVoi    = sort(DoseDis.M240320(indices));
    
    %refGy = round([0.4 0.6 0.8] * max(relevantDose(:)) * 10)/10;
    refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
      
    for runDX = 1:numel(refVol)
       QI.(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));          
    end    
    
    D95.M240320 = QI.D95;        
    HI.M240320 = (DX(2)/DX(98));
    
     % get Dose, dose is sorted to simplify calculations
    relevantDose = DoseDis.WC240320;
    doseInVoi    = sort(DoseDis.WC240320(indices));
    
    %refGy = round([0.4 0.6 0.8] * max(relevantDose(:)) * 10)/10;
    refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
      
    for runDX = 1:numel(refVol)
       QI.(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));          
    end    
    
    D95.WC240320 = QI.D95;        
    HI.WC240320 = (DX(2)/DX(98));
    
        % get Dose, dose is sorted to simplify calculations
    relevantDose = DoseDis.MWC240320;
    doseInVoi    = sort(DoseDis.MWC240320(indices));
    
    %refGy = round([0.4 0.6 0.8] * max(relevantDose(:)) * 10)/10;
    refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
      
    for runDX = 1:numel(refVol)
       QI.(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));          
    end    
    
    D95.MWC240320 = QI.D95;        
    HI.MWC240320 = (DX(2)/DX(98));

figure
plot(weights_Opt)
hold on
plot(weights_WC)
title('Optimized weights for conv. and WC Optimization')

%D95%
disp(D95)
%HI
disp(HI)
    
 %% 4D loop
 % So far the 4D dose distribution was only calculated for a single
 % scenario (
 resultGUI.w = weights_Opt;
 resultGUI.RBExD = DoseDis.Opt240320;
 accDose_Opt = calc4dDose_loop(ct, pln, dij, stf, cst, resultGUI, 'Liver_240320_bf_5mm_new', 'pseudoGTV', 'Leber-GTV', 1 );
 title('static vs moved dose distribution for conventional plan (GTV DVH)')
 
resultGUI.w = weights_WC;
resultGUI.RBExD = DoseDis.WC240320;
accDose_WC = calc4dDose_loop(ct, pln, dij, stf, cst, resultGUI, 'Liver_WC240320_bf_5mm_new', 'pseudoGTV', 'Leber-GTV', 1);
title('static vs moved dose distribution for  worst case plan (GTV DVH)')
 
 