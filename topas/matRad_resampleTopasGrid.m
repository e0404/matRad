function [ctR,cst] = matRad_resampleTopasGrid(ct,cst,pln,stf)
matRad_cfg = MatRad_Config.instance();

if isfield(ct,'modulated') && ct.modulated
    pln.propDoseCalc.useGivenEqDensityCube = 1;
end
matRad_calcDoseInit;

if ~isfield(ct,'resampled')
    
    cubeHUresampled = cell(1,ct.numOfCtScen);
    cubeResampled = cell(1,ct.numOfCtScen);
    for s = 1:ct.numOfCtScen
        cubeHUresampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
        cubeResampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cube{s}, ...
            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
    end
    
    %Allocate temporary resampled CT
    ctR = ct;
    ctR.cube = cell(1);
    ctR.cubeHU = cell(1);
    ctR.numOfCtScen = 1;
    ctR.resolution = dij.doseGrid.resolution;
    ctR.cubeDim = dij.doseGrid.dimensions;
    ctR.x = dij.doseGrid.x;
    ctR.y = dij.doseGrid.y;
    ctR.z = dij.doseGrid.z;
    
    ctR.cubeHU = cubeHUresampled;
    ctR.cube = cubeResampled;

    %Set flag for complete resampling
    ctR.resampled = 1;
else
    ctR = ct;
    matRad_cfg.dispWarning('CT already resampled.');
end

end