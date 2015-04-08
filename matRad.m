% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
%load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
load BOXPHANTOM.mat

% meta information for treatment plan
pln.SAD             = 10000; %[mm]
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0]; % [°]
pln.couchAngles     = [0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = numel(ct.cube);
pln.voxelDimensions = size(ct.cube);
pln.radiationMode   = 'photons'; % either photons / protons / carbon
pln.bioOptimization = false;   % false indicates physical optimization and true indicates biological optimization
pln.numOfFractions  = 30;

% initial visualization
matRad_visCtDose([],cst,pln,ct);

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,0);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
end

%% Dose visualization
doseVis = matRad_mxCalcDose(dij,ones(dij.totalNumOfBixels,1));
matRad_visCtDose(doseVis,cst,pln,ct);

%% change objective function settings if desired
matRad_modCst(cst)

%% inverse planning for imrt
optResult = matRad_inversePlanning(dij,cst,pln);
matRad_visCtDose(optResult,cst,pln,ct);

%% sequencing
if strcmp(pln.radiationMode,'photons')
    Sequencing = matRad_xiaLeafSequencing(optResult.w,stf,7,1);
    seqResult = matRad_mxCalcDose(dij,Sequencing.w);
    matRad_visCtDose(seqResult,cst,pln,ct);
end

%% dvh and conformity index
matRad_calcDVH(optResult,cst)

