function [ct, cst, pln, resultGUI] = matRad_importDicom( files, dicomMetaBool )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad wrapper function to import a predefined set of dicom files into
% matRad's native data formats
% 
% call
%   [ct, cst, pln, resultGUI] = matRad_importDicom( files, dicomMetaBool )
%
% input
%   files:          list of files to be imported (will contain cts and rt
%                   structure set)
%   dicomMetaBool:  (boolean, optional) import complete dicomInfo and
%                   patientName
%
% output
%   ct:        matRad ct struct
%   cst:       matRad cst struct
%   pln:       matRad plan struct
%   resultGUI: matRad result struct holding data for visualization in GUI
%
% References
%   -
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

if ~isdeployed
    addpath('tools');
end

[env, ~] = matRad_getEnvironment();
    
%%
if ~exist('dicomMetaBool','var')
  dicomMetaBool = true;
end

%%
h = waitbar(0,'Please wait...');
%h.WindowStyle = 'Modal';
steps = 2;

%% import ct-cube
waitbar(1 / steps)
resolution.x = files.resx;
resolution.y = files.resy;
resolution.z = files.resz; % [mm] / lps coordinate system
if files.useDoseGrid && isfield(files,'rtdose')
    % get grid from dose cube
    if verLessThan('matlab','9')
        doseInfo = dicominfo(files.rtdose{1,1});
    else
        doseInfo = dicominfo(files.rtdose{1,1},'UseDictionaryVR',true);
    end
    doseGrid{1} = doseInfo.ImagePositionPatient(1) + doseInfo.ImageOrientationPatient(1) * ...
                                                     doseInfo.PixelSpacing(1) * double(0:doseInfo.Columns - 1);
    doseGrid{2} = doseInfo.ImagePositionPatient(2) + doseInfo.ImageOrientationPatient(5) * ...
                                                     doseInfo.PixelSpacing(2) * double(0:doseInfo.Rows - 1);
    doseGrid{3} = doseInfo.ImagePositionPatient(3) + doseInfo.GridFrameOffsetVector(:)';

    % get ct on grid
    ct = matRad_importDicomCt(files.ct, resolution, dicomMetaBool,doseGrid); 

else
    ct = matRad_importDicomCt(files.ct, resolution, dicomMetaBool); 
end

if ~isempty(files.rtss)
    
    %% import structure data
    waitbar(2 / steps)
    structures = matRad_importDicomRtss(files.rtss{1},ct.dicomInfo);
    close(h)

    %% creating structure cube
    h = waitbar(0,'Please wait...');
    %h.WindowStyle = 'Modal';
    steps = numel(structures);
    for i = 1:numel(structures)
        % computations take place here
        waitbar(i / steps)
        fprintf('creating cube for %s volume...\n', structures(i).structName);
        structures(i).indices = matRad_convRtssContours2Indices(structures(i),ct);
    end
    fprintf('finished!\n');
    close(h)

    %% creating cst
    cst = matRad_createCst(structures);

else
    
    cst = matRad_dummyCst(ct);
    
end

%% determine pln parameters
if isfield(files,'rtplan')
    if ~(cellfun(@isempty,files.rtplan(1,:)))
        pln = matRad_importDicomRTPlan(ct, files.rtplan, dicomMetaBool);
    end
end

%% import stf
if isfield(files,'rtplan')
    if ~(cellfun(@isempty,files.rtplan(1,:)))
        if (strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon'))
            %% import steering file
            % pln output because bixelWidth is determined via the stf
            [stf, pln] = matRad_importDicomSteeringParticles(ct, pln, files.rtplan);
        elseif strcmp(pln.radiationMode, 'photons') && isfield(pln.propStf,'collimation')
            % return correct angles in pln 
            [stf, pln] = matRad_importDicomSteeringPhotons(pln);
        else
            warning('No support for DICOM import of steering information for this modality.');
        end
    end
end

%% import dose cube
if isfield(files,'rtdose')
    % check if files.rtdose contains a path and is labeld as RTDose
    % only the first two elements are relevant for loading the rt dose
    if ~(cellfun(@isempty,files.rtdose(1,1:2))) 
        fprintf('loading Dose files \n', structures(i).structName);
        % parse plan in order to scale dose cubes to a fraction based dose
        if exist('pln','var')
            if isfield(pln,'numOfFractions')
                resultGUI = matRad_importDicomRTDose(ct, files.rtdose, pln);
            end
        else
            resultGUI = matRad_importDicomRTDose(ct, files.rtdose);
        end
        if size(resultGUI) == 0
           clear resultGUI;
        end
    end
end

%% put weight also into resultGUI
if exist('stf','var') && exist('resultGUI','var')
    resultGUI.w = [];
    for i = 1:size(stf,2)
        resultGUI.w = [resultGUI.w; [stf(i).ray.weight]'];
    end
end

%% save ct, cst, pln, dose
matRadFileName = [files.ct{1,3} '.mat']; % use default from dicom
[FileName,PathName] = uiputfile('*','Save as...',matRadFileName);
if ischar(FileName)
    % delete unnecessary variables
    switch env
    case 'MATLAB'
        clearvars -except ct cst pln stf resultGUI FileName PathName;
    case 'OCTAVE' 
        clear -x ct cst pln stf resultGUI FileName PathName;
    end
    % save all except FileName and PathName
    save([PathName, FileName], '-regexp', '^(?!(FileName|PathName)$).','-v7.3');
end
