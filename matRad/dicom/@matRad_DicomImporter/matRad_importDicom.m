function matRad_importDicom(obj)
% matRad wrapper function to import a predefined set of dicom files files
% into matRad's native data formats
% 
% In your object, there must be properties that contain: 
%   - list of files to be imported.
% Optional:
%   - а boolean; if you don't want to import complete DICOM information, 
%   set it to false.
%
% Next matRad structures are created in the object and saved in the 
% workspace:
%   - ct, cst, stf, pln, resultGUI.
% *to save them as .mat file you can use matRad_importDicomWidget
%
% call
%   matRad_importDicom(obj)
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();
    
%%
if ~exist('dicomMetaBool','var')
  obj.dicomMetaBool = true;
end

%%
if ~matRad_cfg.disableGUI
    h = waitbar(0,'Please wait...','Color',matRad_cfg.gui.backgroundColor,'DefaultTextColor',matRad_cfg.gui.textColor);
    matRad_applyThemeToWaitbar(h);
else
    h = [];
end
%h.WindowStyle = 'Modal';
steps = 2;

%% import ct-cube
if any(ishandle(h))
    waitbar(1/steps, h)
end
obj.importCT.resolution.x = str2double(obj.importFiles.resx);
obj.importCT.resolution.y = str2double(obj.importFiles.resy);
obj.importCT.resolution.z = str2double(obj.importFiles.resz); % [mm] / lps coordinate system
if obj.importFiles.useImportGrid && isfield(obj.importFiles,'rtdose')
    % get grid from dose cube
    if matRad_cfg.isOctave || verLessThan('matlab','9')
        doseInfo = dicominfo(obj.importFiles.rtdose{1,1});
    else
        doseInfo = dicominfo(obj.importFiles.rtdose{1,1},'UseDictionaryVR',true);
    end
    obj.ImportGrid{1} = doseInfo.ImagePositionPatient(1) + doseInfo.ImageOrientationPatient(1) * ...
                                                     doseInfo.PixelSpacing(1) * double(0:doseInfo.Columns - 1);
    obj.ImportGrid{2} = doseInfo.ImagePositionPatient(2) + doseInfo.ImageOrientationPatient(5) * ...
                                                     doseInfo.PixelSpacing(2) * double(0:doseInfo.Rows - 1);
    obj.ImportGrid{3} = doseInfo.ImagePositionPatient(3) + doseInfo.GridFrameOffsetVector(:)';

    % get ct on grid
    obj = matRad_importDicomCt(obj); 

else
    obj = matRad_importDicomCt(obj); 
end

if ~isempty(obj.importFiles.rtss)
    
    %% import structure data
    if any(ishandle(h))
        waitbar(2/steps, h)
    end
    obj = matRad_importDicomRtss(obj);
    if any(ishandle(h))
        close(h)
    end

    %% creating structure cube
    if ~matRad_cfg.disableGUI
        h = waitbar(0,'Please wait...','Color',matRad_cfg.gui.backgroundColor,'DefaultTextColor',matRad_cfg.gui.textColor);
        matRad_applyThemeToWaitbar(h);
    else
        h = [];
    end
    %h.WindowStyle = 'Modal';
    steps = numel(obj.importRtss.structures);


    % The x- & y-direction in lps-coordinates are specified in:
    % ImageOrientationPatient

obj.importRtss.xDir = obj.ct.dicomInfo.ImageOrientationPatient(1:3); % lps: [1;0;0]
obj.importRtss.yDir = obj.ct.dicomInfo.ImageOrientationPatient(4:6); % lps: [0;1;0]


if ~(obj.importRtss.xDir(1) == 1 && obj.importRtss.xDir(2) == 0 && obj.importRtss.xDir(3) == 0)
     matRad_cfg.dispInfo('\nNonstandard image orientation: tring to Mirror RTSS x-direction...')
end

if ~(obj.importRtss.yDir(1) == 0 && obj.importRtss.yDir(2) == 1 && obj.importRtss.yDir(3) == 0)
    matRad_cfg.dispInfo('\nNonstandard image orientation: trying to Mirror RTSS y direction...')
end
    for i = 1:numel(obj.importRtss.structures)
        % computations take place here
        if any(ishandle(h))
            waitbar(1/steps, h)
        end
        matRad_cfg.dispInfo('creating cube for %s volume... ', obj.importRtss.structures(i).structName);
        try
            obj.importRtss.structures(i).indices = matRad_convRtssContours2Indices(obj.importRtss.structures(i),obj.ct);
            matRad_cfg.dispInfo('\n');
        catch ME
            warning('matRad:dicomImport','could not be imported: %s',ME.message);
            obj.importRtss.structures(i).indices = [];
        end      
    end
    matRad_cfg.dispInfo('finished!\n');
    close(h)

    %% creating cst
    obj = matRad_createCst(obj);

else
    
    obj = matRad_dummyCst(obj);
    
end

%% determine pln parameters
if ~isempty(obj.importFiles.rtplan)
    if ~(cellfun(@isempty,obj.importFiles.rtplan(1,:)))
        obj = matRad_importDicomRTPlan(obj);
    end
else
    obj.pln = struct([]);
end

%% import stf
if ~isempty(obj.importFiles.rtplan)
    if ~(cellfun(@isempty,obj.importFiles.rtplan(1,:)))
        if (strcmp(obj.pln.radiationMode,'protons') || strcmp(obj.pln.radiationMode,'carbon'))
            %% import steering file
            % pln output because bixelWidth is determined via the stf
            obj = matRad_importDicomSteeringParticles(obj);
        elseif strcmp(obj.pln.radiationMode, 'photons') && isfield(obj.pln.propStf,'collimation')
            % return correct angles in pln 
            obj = matRad_importDicomSteeringPhotons(obj);
        else
            warning('No support for DICOM import of steering information for this modality.');
        end
    end
else
    obj.stf = struct([]);
end

%% import dose cube
if ~isempty(obj.importFiles.rtdose) 
    % check if obj.importFiles.rtdose contains a path and is labeld as RTDose
    % only the first two elements are relevant for loading the rt dose
    if ~(cellfun(@isempty,obj.importFiles.rtdose(1,1))) 
        matRad_cfg.dispInfo('loading dose files...\n');
        % parse plan in order to scale dose cubes to a fraction based dose
        obj = matRad_importDicomRTDose(obj);        
        if size(obj.resultGUI) == 0
           obj.resultGUI = struct([]);
        end
    end
    matRad_cfg.dispInfo('finished!\n');
else
    obj.resultGUI = struct([]);
    matRad_cfg.dispInfo('There are no dose files!\n');
end

%% put weight also into resultGUI
if ~isempty(obj.stf) && ~isempty(obj.resultGUI)
    obj.resultGUI.w = [];
    for i = 1:size(obj.stf,2)
        obj.resultGUI.w = [obj.resultGUI.w; [obj.stf(i).ray.weight]'];
    end
end
if any(ishandle(h))
        close(h)
end
%% put ct, cst, pln, stf, resultGUI to the workspace
ct = obj.ct;
cst = obj.cst;
pln = obj.pln;
stf = obj.stf;
resultGUI = obj.resultGUI;

assignin('base', 'ct', ct);
assignin('base', 'cst', cst);

if ~isempty(obj.pln)
    assignin('base', 'pln', pln);
end 

if ~isempty(obj.stf)
    assignin('base', 'stf', stf);
end

if ~isempty(obj.resultGUI)
assignin('base', 'resultGUI', resultGUI);
end

end


