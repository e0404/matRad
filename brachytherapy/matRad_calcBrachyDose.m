function dij = matRadBrachy_calcDose(ct,stf,pln,cst)
% Dose calculation
%%Configure
matRad_cfg =  MatRad_Config.instance();
matRad_calcBrachyDoseInit;


%% get dose points
[XGrid,YGrid,ZGrid] = meshgrid(dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
dosePoints.x = reshape(XGrid,1,[]);
dosePoints.y = reshape(YGrid,1,[]);
dosePoints.z = reshape(ZGrid,1,[]);

%% get seed points
seedPoints.x = reshape(stf.seedPosX,1,[]);
seedPoints.y = reshape(stf.seedPosY,1,[]);
seedPoints.z = reshape(stf.seedPosZ,1,[]);

%% get seed dosepoint distance matrix
[DistanceMatrix,DistanceVector] = getDistanceMatrix(seedPoints,dosePoints);

%% seed dosepoint angle matrix
if strcmp(pln.propDoseCalc.TG43approximation,'2D')
    templateNormal = normalize(pln.propStf.orientation.Zdir,'norm');
    [ThetaMatrix,ThetaVector] = getThetaMatrix(templateNormal,DistanceMatrix);
end


%% Calculate Dose Rate matrix
isotope = Source(dij.basedata.data); %create instance of dose calculation class

switch pln.propDoseCalc.TG43approximation
    case '1D'
        Dose = isotope.getDoseRate1D(DistanceVector);
    case '2D'
        Dose = isotope.getDoseRate2D(DistanceVector,ThetaVector);
end
        
dij.physicalDose = {sparse(reshape(Dose,size(DistanceMatrix.norm)))};
end



