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
load LIVER.mat
%load BOXPHANTOM.mat

% meta information for treatment plan
pln.SAD             = 1000; %[mm]
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0:72:359]; % [°]
pln.couchAngles     = [0 0 0 0 0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = numel(ct.cube);
pln.voxelDimensions = size(ct.cube);
pln.radiationMode   = 'photons'; % either photons / protons / carbon
pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 1;

%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,0);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
end

%% inverse planning for imrt
resultGUI = matRad_fluenceOptimization(dij,cst,pln,0);

%% sequencing
if strcmp(pln.radiationMode,'photons')
    %Sequencing = matRad_xiaLeafSequencing(resultGUI.w,stf,7,1);
    Sequencing = matRad_engelLeafSequencing(resultGUI.w,stf,7);
    resultGUI = matRad_mxCalcDose(dij,Sequencing.w,cst);
end

%% start gui for visualization of result
matRadGUI

%% dvh
matRad_calcDVH(resultGUI,cst)
