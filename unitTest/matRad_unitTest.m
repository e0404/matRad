% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad unit test script
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
clc, clear, close all

%% meta information for treatment plan - new test scenarios can be created here

pln.SAD             = 1000; %[mm]
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0]; % [°]
pln.couchAngles     = [0]; % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.radiationMode   = 'protons'; % either photons / protons / carbon
pln.bioOptimization = 'none'; % none: physical optimization; effect: effect-based optimization; RBExD: optimization of RBE-weighted dose
pln.numOfFractions  = 20;
pln.runSequencing   = true; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = true; % 1/true: run DAO, 0/false: don't / will be ignored for particles

FlagSuccess = matRad_unitTestGenData(pln,'BOXPHANTOM');

%% run tests
addpath(fullfile(pwd,'..'));
SaveResultToDisk = true;
% test first scenario
matRad_unitTestRun('unitTest_result_X_001_BOXPHANTOM',SaveResultToDisk);

% test second scenario
matRad_unitTestRun('unitTest_result_P_001_BOXPHANTOM',SaveResultToDisk); 



