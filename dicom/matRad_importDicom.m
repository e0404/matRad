function [ct,cst,pln,stf,resultGUI] = matRad_importDicom( files, dicomMetaBool )
% matRad wrapper function to import a predefined set of dicom files 
% into matRad's native data formats
% 
% call
% [ct, cst, pln, stf, resultGUI] = matRad_importDicom( files )  
% [ct, cst, pln, stf, resultGUI] = matRad_importDicom( files, dicomMetaBool )
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
%   stf:       matRad stf struct
%   resultGUI: matRad result struct holding data for visualization in GUI
%
% References
%   -
%
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
matRad_cfg = MatRad_Config.instance();
    
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
    if matRad_cfg.isOctave || verLessThan('matlab','9')
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


    % The x- & y-direction in lps-coordinates are specified in:
    % ImageOrientationPatient

xDir = ct.dicomInfo.ImageOrientationPatient(1:3); % lps: [1;0;0]
yDir = ct.dicomInfo.ImageOrientationPatient(4:6); % lps: [0;1;0]


if ~(xDir(1) == 1 && xDir(2) == 0 && xDir(3) == 0)
     matRad_cfg.dispInfo('\nNonstandard image orientation: tring to Mirror RTSS x-direction...')
end

if ~(yDir(1) == 0 && yDir(2) == 1 && yDir(3) == 0)
    matRad_cfg.dispInfo('\nNonstandard image orientation: trying to Mirror RTSS y direction...')
end
    for i = 1:numel(structures)
        % computations take place here
        waitbar(i / steps)
        fprintf('creating cube for %s volume... ', structures(i).structName);
        try
            structures(i).indices = matRad_convRtssContours2Indices(structures(i),ct);
            fprintf('\n');
        catch ME
            warning('matRad:dicomImport','could not be imported: %s',ME.message);
            structures(i).indices = [];
        end      
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
else
    pln = struct([]);
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
else
    stf = struct([]);
end

%% import dose cube
if isfield(files,'rtdose')
    % check if files.rtdose contains a path and is labeld as RTDose
    % only the first two elements are relevant for loading the rt dose
    if ~(cellfun(@isempty,files.rtdose(1,1:2))) 
        fprintf('loading Dose files \n', structures(i).structName);
        % parse plan in order to scale dose cubes to a fraction based dose
        if exist('pln','var') && ~isempty(pln) && isfield(pln,'numOfFractions')
            resultGUI = matRad_importDicomRTDose(ct, files.rtdose, pln);
        else
            resultGUI = matRad_importDicomRTDose(ct, files.rtdose);
        end
        if size(resultGUI) == 0
           resultGUI = struct([]);
        end
    end
else
    resultGUI = struct([]);
end

%% put weight also into resultGUI
if ~isempty(stf) && ~isempty(resultGUI)
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
    if matRad_cfg.isMatlab
        varNames=who;
        for var=1:length(varNames)
            cleanVar=eval([varNames{var}]);
            if isempty(cleanVar)
                eval(['clear ' varNames{var} ';'])
            end
        end
        clearvars -except ct cst pln stf resultGUI FileName PathName;
        save([PathName, FileName], '-regexp', '^(?!(FileName|PathName)$).','-v7');
    elseif matRad_cfg.isOctave
        varNames=who;
        for var=1:length(varNames)
            cleanVar=eval([varNames{var}]);
            if isempty(cleanVar)
                eval(['clear ' varNames{var} ';'])
            end
        end
        clear -x ct cst pln stf resultGUI FileName PathName;
        save([PathName, FileName],'-v6');
    else
    end
    % save all except FileName and PathName
    
end
