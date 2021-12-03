function dij = matRad_calcBrachyDose(ct,stf,pln,cst)
% matRad_calcBrachyDose calculates dose influence matrix according to the
% AAPM update Rivard et al. 2004
%
% call
%   dij = matRad_calcBrachyDose(ct,stf,pln,cst)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%               (positions and constraints of patient structures)
%   pln:        matRad plan meta information struct
%   stf:        struct containing geometric information
%
% output
%   dij:        stuct containing dose influence information
%
% References: 
%   [1] https://doi.org/10.1118/1.1646040 - TG43 Update
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Configure
matRad_cfg =  MatRad_Config.instance();
matRad_calcBrachyDoseInit;

% initialize waitbar (always indented to seperate from important code)
figureWait = waitbar...
    (0,'calculating dose inlfluence matrix for brachytherapy...');
matRad_cfg.dispInfo('Starting  brachytherapy dose calculation...\n');
startTime = tic;

%% get dose points and seedpoints
% "dosePoints" and "seedPoints" are both structs with fields x,y,z:
% each contains a 1D row vector of position components [mm]

seedPoints.x = single(stf.seedPoints.x);
seedPoints.y = single(stf.seedPoints.y);
seedPoints.z = single(stf.seedPoints.z);

[XGrid,YGrid,ZGrid] = meshgrid(dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
dosePoints.x = single(reshape(XGrid,1,[]));
dosePoints.y = single(reshape(YGrid,1,[]));
dosePoints.z = single(reshape(ZGrid,1,[]));

matRad_cfg.dispInfo('\t computing distance transform... ');


%% get seed dosepoint distance matrix
% [seedPoint x dosePoint] matrix with relative distance as entries
% detailed documentation in function
DistanceMatrix = matRad_getDistanceMatrix(seedPoints,dosePoints);

% ignore all distances > Cutoff for the following calculations to save time
Ignore = DistanceMatrix.dist > pln.propDoseCalc.DistanceCutoff;
calcDistanceMatrix.x = DistanceMatrix.x(~Ignore);
calcDistanceMatrix.y = DistanceMatrix.y(~Ignore);
calcDistanceMatrix.z = DistanceMatrix.z(~Ignore);
calcDistanceMatrix.dist = DistanceMatrix.dist(~Ignore);

% now all fields of calcDistanceMatrix are n x 1 arrays!

% update waitbar
waitbar(0.125);
matRad_cfg.dispInfo('done in %f s!\n',toc(startTime));

%% seed dosepoint angle matrix
% [seedPoint x dosePoint] matrix with relative theta angle as entries
% detailed documentation in function
% only call for 2D formalism
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'TG43approximation')
    pln.propDoseCalc.TG43approximation = '2D';
end

if strcmp(pln.propDoseCalc.TG43approximation,'2D')
    matRad_cfg.dispInfo('\t computing angle for TG43-2D... ');
    tmpTimer = tic;
    [ThetaMatrix,~] = matRad_getThetaMatrix(pln.propStf.template.normal,calcDistanceMatrix);
    matRad_cfg.dispInfo('done in %f s!\n',toc(tmpTimer));
end

% update waitbar
waitbar(0.25);
%% Calculate Dose Rate matrix
% Calculation according to [1]

matRad_cfg.dispInfo('\t computing dose-rate for TG43-%s... ',pln.propDoseCalc.TG43approximation);
tmpTimer = tic;
DoseRate = zeros(length(dosePoints.x),length(seedPoints.x));
switch pln.propDoseCalc.TG43approximation
    case '1D'        
        DoseRate(~Ignore) = ...
        matRad_getDoseRate1D_poly(machine,calcDistanceMatrix.dist);
    case '2D'
        DoseRate(~Ignore) = ...
        matRad_getDoseRate2D_poly(machine,calcDistanceMatrix.dist,ThetaMatrix);
    otherwise
        matRad_cfg.dispError('TG43 Approximation ''%s'' not known!',pln.propDoseCalc.TG43approximation);
end
matRad_cfg.dispInfo('done in %f s!\n',toc(tmpTimer));

dij.physicalDose = {DoseRate};

% update waitbar, delete waitbar
waitbar(1);
matRad_cfg.dispInfo('Brachytherapy dose calculation finished in %f s!\n',toc(startTime));
delete(figureWait);


end



