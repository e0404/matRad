function dij = matRad_calcBrachyDose(ct,stf,pln,cst)
% matRad_calcBrachyDose calculates dose influence matrix according to the
% AAPM update Rivard et al. 2004
%
% call
%   dij = matRad_calcBrachyDose(ct,stf,pln,cst)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct (positions and constraints of patient structures)
%   pln:        matRad plan meta information struct
%   stf:        struct containing geometric information about dose and seed points
%
% output
%   dij:        stuct containing dose influence information

tic;

%%Configure
matRad_cfg =  MatRad_Config.instance();
matRad_calcBrachyDoseInit;



%% get dose points and seedpoints
% "dosePoints" and "seedPoints" are both structs with fields x,y,z:
% each contains a 1D row vector of position components 
[XGrid,YGrid,ZGrid] = meshgrid(dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
dosePoints.x = reshape(XGrid,1,[]);
dosePoints.y = reshape(YGrid,1,[]);
dosePoints.z = reshape(ZGrid,1,[]);

seedPoints = stf.seedPoints;
%% get seed dosepoint distance matrix
% [seedPoint x dosePoint] matrix with relative distance as entries
% detailed documentation in function
[DistanceMatrix,DistanceVector] = getDistanceMatrix(seedPoints,dosePoints);
toc;
%% seed dosepoint angle matrix
% [seedPoint x dosePoint] matrix with relative theta angle as entries
% detailed documentation in function
% only call for 2D formalism
if strcmp(pln.propDoseCalc.TG43approximation,'2D')
    Zdir = pln.propStf.shiftRotMtx(1:3,3);
    templateNormal = normalize(Zdir,'norm');
    [ThetaMatrix,ThetaVector] = getThetaMatrix(templateNormal,DistanceMatrix);
end

toc;
%% Calculate Dose Rate matrix
% Calculation according to  Rivard et al. (2004): AAPM TG-43 update
% implemented by Dr. Christian Guthier

%create instance of dose calculation class input: machine basedata. 
% For required fields and their format see matRad/basedata/brachy_LDR or
% class def
switch pln.propDoseCalc.TG43approximation
    case '1D'
        DoseRate = matRad_getDoseRate1D_poly(machine,DistanceMatrix.dist,ThetaMatrix);
    case '2D'
        DoseRate = matRad_getDoseRate2D_poly(machine,DistanceMatrix.dist,ThetaMatrix);
end
toc;
dij.physicalDose = {DoseRate};
toc;
end



