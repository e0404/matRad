numLevels = 1;
numApertures = 3:2:pln.maxNumApertures;

CIFMO = zeros(numel(numLevels),numel(numApertures));
HIFMO = zeros(numel(numLevels),numel(numApertures));

planModulationSeq = zeros(numel(numLevels),numel(numApertures));
planMUSeq = zeros(numel(numLevels),numel(numApertures));
planAreaSeq = zeros(numel(numLevels),numel(numApertures));
CISeq = zeros(numel(numLevels),numel(numApertures));
HISeq = zeros(numel(numLevels),numel(numApertures));

planModulationSeq2 = zeros(numel(numLevels),numel(numApertures));
planMUSeq2 = zeros(numel(numLevels),numel(numApertures));
planAreaSeq2 = zeros(numel(numLevels),numel(numApertures));
planTimeSeq2 = zeros(numel(numLevels),numel(numApertures));
CISeq2 = zeros(numel(numLevels),numel(numApertures));
HISeq2 = zeros(numel(numLevels),numel(numApertures));

planModulation = zeros(numel(numLevels),numel(numApertures));
planMU = zeros(numel(numLevels),numel(numApertures));
planArea = zeros(numel(numLevels),numel(numApertures));
planTimeBeforeOpt = zeros(numel(numLevels),numel(numApertures));
planTime = zeros(numel(numLevels),numel(numApertures));
CI = zeros(numel(numLevels),numel(numApertures));
HI = zeros(numel(numLevels),numel(numApertures));

i = 1;
j = 1;

for numAp = numApertures
    %change number of apertures, changes initialized and optimized angles
    pln.numApertures = numAp;
    pln.VMAT = 1;
    pln = matRad_VMATGantryAngles(pln,'new');
    fprintf('\nNumber of apertures = %d\n',numAp);
    
    %update steering file
    stf = matRad_generateStf(ct,cst,pln);
    
    %dij doesn't need updating since the gantry angles are always the same
    
    %do fluence optimization on new initialized angles
    resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf);
    
    %determine quality indicators
    resultGUI = matRad_calcQualityIndicators(resultGUI,cst,pln);
    
    %enter QI
    %CIFMO(:,j) = resultGUI.QI(3).CI_1_67Gy;
    %HIFMO(:,j) = resultGUI.QI(3).HI_1_67Gy;
    
    %H&N
    CIFMO(:,j) = resultGUI.QI(19).CI_2_33Gy;
    HIFMO(:,j) = resultGUI.QI(19).HI_2_33Gy;

    %Prostate
    %CIFMO(:,j) = resultGUI.QI(6).CI_2_27Gy;
    %HIFMO(:,j) = resultGUI.QI(6).HI_2_27Gy;
    
    %do leaf sequencing, first keeping all apertures
    pln.VMAT = 0;
    resultGUI = matRad_svenssonLeafSequencing(resultGUI,stf,dij,pln,numAp,0);
    
    %determine plan metrics
    resultGUI = matRad_calcDeliveryMetrics(resultGUI,pln,stf);
    resultGUI = matRad_calcQualityIndicators(resultGUI,cst,pln);
    
    %enter metrics
    planModulationSeq(i,j) = resultGUI.apertureInfo.planModulation;
    planMUSeq(i,j) = resultGUI.apertureInfo.planMU;
    planAreaSeq(i,j) = resultGUI.apertureInfo.planArea;
    %CISeq(i,j) = resultGUI.QI(3).CI_1_67Gy;
    %HISeq(i,j) = resultGUI.QI(3).HI_1_67Gy;
    
    %H&N
    CISeq(i,j) = resultGUI.QI(19).CI_2_33Gy;
    HISeq(i,j) = resultGUI.QI(19).HI_2_33Gy;
    
    %Prostate
    %CISeq(i,j) = resultGUI.QI(6).CI_2_27Gy;
    %HISeq(i,j) = resultGUI.QI(6).HI_2_27Gy;
    
    
    
    %now do leaf sequencing, but preparing for DAO (spread apertures)
    pln.VMAT = 1;
    resultGUI = matRad_svenssonLeafSequencing(resultGUI,stf,dij,pln,numAp,0);
    
    %determine plan metrics
    resultGUI = matRad_calcDeliveryMetrics(resultGUI,pln,stf);
    resultGUI = matRad_calcQualityIndicators(resultGUI,cst,pln);
    
    %enter metrics
    planModulationSeq2(i,j) = resultGUI.apertureInfo.planModulation;
    planMUSeq2(i,j) = resultGUI.apertureInfo.planMU;
    planAreaSeq2(i,j) = resultGUI.apertureInfo.planArea;
    planTimeSeq2(i,j) = resultGUI.apertureInfo.time;
    %CISeq2(i,j) = resultGUI.QI(3).CI_1_67Gy;
    %HISeq2(i,j) = resultGUI.QI(3).HI_1_67Gy;
    
    %H&N
    CISeq2(i,j) = resultGUI.QI(19).CI_2_33Gy;
    HISeq2(i,j) = resultGUI.QI(19).HI_2_33Gy;
    
    %Prostate
    %CISeq2(i,j) = resultGUI.QI(6).CI_2_27Gy;
    %HISeq2(i,j) = resultGUI.QI(6).HI_2_27Gy;
    
    
    
    %do DAO
    resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,1);
    
    resultGUI = matRad_calcDeliveryMetrics(resultGUI,pln,stf);
    planTimeBeforeOpt(i,j) = resultGUI.apertureInfo.time;
    
    %optimize delivery
    resultGUI = matRad_optDelivery(resultGUI,pln,1);
    
    %determine plan metrics
    resultGUI = matRad_calcDeliveryMetrics(resultGUI,pln,stf);
    resultGUI = matRad_calcQualityIndicators(resultGUI,cst,pln);
    
    %enter metrics
    planModulation(i,j) = resultGUI.apertureInfo.planModulation;
    planMU(i,j) = resultGUI.apertureInfo.planMU;
    planArea(i,j) = resultGUI.apertureInfo.planArea;
    planTime(i,j) = resultGUI.apertureInfo.time;
    %CI(i,j) = resultGUI.QI(3).CI_1_67Gy;
    %HI(i,j) = resultGUI.QI(3).HI_1_67Gy;
    
    %H&N
    CI(i,j) = resultGUI.QI(19).CI_2_33Gy;
    HI(i,j) = resultGUI.QI(19).HI_2_33Gy;
    
    %Prostate
    %CI(i,j) = resultGUI.QI(6).CI_2_27Gy;
    %HI(i,j) = resultGUI.QI(6).HI_2_27Gy;
    
    fname = sprintf('H&N Svensson Apertures %d',numAp);
    save(fname,'resultGUI');
    
    close all
    
    j = j+1;
end

save('Results H&N Svensson','plan*','num*','CI*','HI*')


