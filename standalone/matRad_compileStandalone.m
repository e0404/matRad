% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to compile an executable matRad version with the matlab
% application compiler.
%
% There will be no splash screen included in the executable if following 
% the instructions below.
% In order to compile an executable with splash screen call deploytool and
% add the required files manually.
%
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

    
%% start the compiling process from the matRad root directory
dirMatRad = [pwd filesep];                                    % set the matRad root directory
dirOut = [pwd  filesep 'standalone' filesep 'build' filesep]; % specify the output folder

if exist(dirOut)~= 7
   mkdir(dirOut);
end

addpath('dicomImport');         % include the DICOM Import files in the dependency analysis

mcc('-C','-o','matRad','-W','WinMain:matRad','-T','link:exe','-d',dirOut,...      
    '-R','-logfile,consoleOutput.log','-v',[dirMatRad 'matRadGUI.m'],...
    '-a',[dirMatRad 'dicomImport' filesep 'DKFZ_Logo.png'],...
    '-a',[dirMatRad 'dicomImport' filesep 'matrad_logo.png'],...
    '-a',[dirMatRad 'dicomImport' filesep 'hlutDefault.mat'],...
    '-a',[dirMatRad 'photons_Generic.mat'],...
    '-a',[dirMatRad 'protons_Generic.mat'],...
    '-a',[dirMatRad 'carbon_Generic.mat'])
