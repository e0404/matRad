DmaxVec = 5:5:25;
wCoreVec = 100:100:1000;

CIFMO = zeros(numel(DmaxVec),numel(wCoreVec));
D10FMO = zeros(numel(DmaxVec),numel(wCoreVec));

CISeq = zeros(numel(DmaxVec),numel(wCoreVec));
D10Seq = zeros(numel(DmaxVec),numel(wCoreVec));

CISeq2 = zeros(numel(DmaxVec),numel(wCoreVec));
D10Seq2 = zeros(numel(DmaxVec),numel(wCoreVec));

planTimeBeforeOpt = zeros(numel(DmaxVec),numel(wCoreVec));
planTime = zeros(numel(DmaxVec),numel(wCoreVec));
CI = zeros(numel(DmaxVec),numel(wCoreVec));
D10 = zeros(numel(DmaxVec),numel(wCoreVec));

i = 1;
j = 1;

for wCore = wCoreVec
    for Dmax = DmaxVec
        cst{2,6}.penalty = wCore;
        cst{2,6}.dose = Dmax;
        pln.VMAT = 1;
        pln = matRad_VMATGantryAngles(pln,'new');
        fprintf('\nCore weight = %d, max Dose = %d Gy\n',wCore,Dmax);
        
        %dij doesn't need updating since the gantry angles are always the same
        
        %do fluence optimization on new initialized angles
        resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf,1);
        
        %determine quality indicators
        resultGUI = matRad_calcQualityIndicators(resultGUI,cst,pln);
        
        %enter QI
        CIFMO(i,j) = resultGUI.QI(3).CI_2Gy;
        D10FMO(i,j) = resultGUI.QI(2).D10;
        
        
        fname = sprintf('TG119 D10 wCore %d, Dmax %d',wCore,Dmax);
        save(fname,'resultGUI');
        
        i = i+1;
    end
    i = 1;
    j = j+1;
end

save('Results TG119 D10','DmaxVec','wCoreVec','CI*','D10*')

