function [ctR,cst,stf] = matRad_resampleCTtoGrid(ct,cst,pln,stf)
% function to resample the ct grid for example for faster MC computation
%
% call
%   [ctR,cst,stf] = matRad_resampleGrid(ct,cst,stf)
%
% input
%   ct:             Path to folder where TOPAS files are in (as string)
%   cst:            matRad segmentation struct
%   pln:            matRad plan struct
%   stf:            matRad steering struct
%
% output
%   ctR:            resampled CT
%   cst:            updated ct struct (due to calcDoseInit)
%   stf:            updated stf struct (due to calcDoseInit)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class

% Load calcDoseInit as usual
% Note of warning, even though the pln is marked as unused, it is needed for calcDoseInit!
matRad_calcDoseInit;

% Check if CT has already been resampled
if ~isfield(ct,'resampled')
    % Allpcate resampled cubes
    cubeHUresampled = cell(1,ct.numOfCtScen);
    cubeResampled = cell(1,ct.numOfCtScen);

    % Perform resampling to dose grid
    for s = 1:ct.numOfCtScen
        cubeHUresampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{s}, ...
            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
        cubeResampled{s} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cube{s}, ...
            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');
    end

    % Allocate temporary resampled CT
    ctR = ct;
    ctR.cube = cell(1);
    ctR.cubeHU = cell(1);

    % Set CT resolution to doseGrid resolution
    ctR.resolution = dij.doseGrid.resolution;
    ctR.cubeDim = dij.doseGrid.dimensions;
    ctR.x = dij.doseGrid.x;
    ctR.y = dij.doseGrid.y;
    ctR.z = dij.doseGrid.z;

    % Write resampled cubes
    ctR.cubeHU = cubeHUresampled;
    ctR.cube = cubeResampled;

    % Set flag for complete resampling
    ctR.resampled = 1;
    ctR.ctGrid = dij.doseGrid;

    % Save original grid
    ctR.originalGrid = dij.ctGrid;
else
    ctR = ct;
    matRad_cfg.dispWarning('CT already resampled.');
end
end
