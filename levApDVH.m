% script to calculate DVH metrics as a function of number of levels and
% number of kept apertures

sprintf('Script will calculate DVH metrics for sequenced plan, varying the number of levels and number of apertures to keep');

str = input('Select min number of levels to use: ','s');
minNumLev = str2double(str);
str = input('Select max number of levels to use: ','s');
maxNumLev = str2double(str);

str = input('Select min number of apertures to keep: ','s');
minApKeep = str2double(str);
str = input('Select max number of apertures to keep: ','s');
maxApKeep = str2double(str);

%recalculate dose: optimal fluence
resultGUI.physicalDose = reshape(dij.physicalDose{1} * resultGUI.wUnsequenced,dij.dimensions);
resultGUI = matRad_calcQualityIndicators(resultGUI,cst,pln);
refMeanDose = resultGUI.QI(3).mean;
refStd = resultGUI.QI(3).std;
refD98 = resultGUI.QI(3).D98;
refD95 = resultGUI.QI(3).D95;
refCI = resultGUI.QI(3).CI_1_67Gy;
refHI = resultGUI.QI(3).HI_1_67Gy;

meanDose = zeros(maxNumLev-minNumLev+1,maxApKeep-minApKeep+1);
std = meanDose;
D98 = meanDose;
D95 = meanDose;
CI = meanDose;
HI = meanDose;
ap = meanDose;
numberOfLevels = meanDose;
numberOfApertures = meanDose;

i = 1;
optCI = 0;
optHI = 100;
for numLev = minNumLev:maxNumLev
    j = 1;
    for apKeep = minApKeep:maxApKeep
        
        resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,numLev,0,0,apKeep);
        resultGUI = matRad_calcQualityIndicators(resultGUI,cst,pln);
        
        meanDose(i,j) = resultGUI.QI(3).mean;
        std(i,j) = resultGUI.QI(3).std;
        D98(i,j) = resultGUI.QI(3).D98;
        D95(i,j) = resultGUI.QI(3).D95;
        CI(i,j) = resultGUI.QI(3).CI_1_67Gy;
        HI(i,j) = resultGUI.QI(3).HI_1_67Gy;
        numberOfLevels(i,j) = numLev;
        numberOfApertures(i,j) = apKeep;
        
        if CI(i,j) >= optCI && apKeep ~= 0
            optCI = CI(i,j);
            optCINumLev = numLev;
            optCIApKeep = apKeep;
        end
        if HI(i,j) <= optHI && apKeep ~= 0
            optHI = HI(i,j);
            optHINumLev = numLev;
            optHIApKeep = apKeep;
        end
        
        j = j+1;
    end
    i = i+1;
end