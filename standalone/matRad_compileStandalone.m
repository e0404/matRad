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
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%% start the compiling process from the matRad root directory
dirMatRad = [pwd filesep];                                    % set the matRad root directory
dirOut = [pwd  filesep 'standalone' filesep 'build' filesep]; % specify the output folder

if exist(dirOut)~= 7
   mkdir(dirOut);
end

addpath('dicomImport');         % include the DICOM Import files in the dependency analysis
addpath('optimization');         % include the DICOM Import files in the dependency analysis

mcc('-C','-o','matRad','-W','WinMain:matRad','-T','link:exe','-d',dirOut,...      
    '-R','-logfile,matRad_consoleOutput.log','-v',[dirMatRad 'matRadGUI.m'],...
    '-a',[dirMatRad 'dicomImport' filesep 'DKFZ_Logo.png'],...
    '-a',[dirMatRad 'dicomImport' filesep 'matrad_logo.png'],...
    '-a',[dirMatRad 'dicomImport' filesep 'hlutLibrary\matRad_default.hlut'],...
    '-a',[dirMatRad 'photons_Generic.mat'],...
    '-a',[dirMatRad 'protons_Generic.mat'],...
    '-a',[dirMatRad 'carbon_Generic.mat'])


