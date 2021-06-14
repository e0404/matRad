function dij = matRadBrachy_calcDose(ct,stf,pln,cst)
% Dose calculation
%%Configure
switch(pln.machine)
    case 'Generic'
        pln.machine = 'Amersham6711';
matRad_cfg =  MatRad_Config.instance();
matRad_calcBrachyDoseInit;
end

%% get dose points
[XGrid,YGrid,ZGrid] = meshgrid(dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
XdosePoints = reshape(XGrid,1,[]);
YdosePoints = reshape(YGrid,1,[]);
ZdosePoints = reshape(ZGrid,1,[]);

%% get seed points
XseedPoints = reshape(stf.seedPosX,1,[]);
YseedPoints = reshape(stf.seedPosY,1,[]);
ZseedPoints = reshape(stf.seedPosZ,1,[]);

%% get seed dosepoint distance matrix
% distance vector matrix
Xdiff = XdosePoints'*ones(1,length(XseedPoints)) - ones(length(XdosePoints),1)*XseedPoints;
Ydiff = YdosePoints'*ones(1,length(YseedPoints)) - ones(length(YdosePoints),1)*YseedPoints;
Zdiff = ZdosePoints'*ones(1,length(ZseedPoints)) - ones(length(ZdosePoints),1)*ZseedPoints;

DistanceMatrix = sqrt(Xdiff.^2+Ydiff.^2+Zdiff.^2);
%% seed dosepoint angle matrix
% only for 2D...

%% Calculate Dose Rate matrix
physicalDose = getPointDose1D(dij.basedata,DistanceMatrix,pln.propDoseCalc.durationImplanted);
 dij.physicalDose = {sparse(physicalDose)};
end



