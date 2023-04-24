addpath(fileparts(pwd));
matRad_rc
load referenceDoseCalculation


%% 1D

    machine = refDos.TG43_1D.basedata;
    r = refDos.coords.r; 

    doseCal = matRad_getDoseRate1D_poly(machine,r);
    doseRef = refDos.TG43_1D.fullDose;
    
    [gammaPassRate,~] = ...
    matRad_gammaPassRate2D(doseRef,doseCal, grids, 3 , 3 , 0.1);

%% 2D
    machine = refDos.TG43_2D.basedata;
    r = refDos.coords.r; 
    theta = refDos.coords.theta; 

    doseCal = matRad_getDoseRate2D_poly(machine,r,theta);
    doseRef = refDos.TG43_2D.fullDose;
    
    [gammaPassRate,~] = ...
    matRad_gammaPassRate2D(doseRef,doseCal, grids, 3 , 3 , 0.1);

