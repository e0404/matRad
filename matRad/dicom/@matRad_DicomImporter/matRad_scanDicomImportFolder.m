function obj = matRad_scanDicomImportFolder(obj)
% matRad function to scan a folder for dicom data
% 
% call
%   obj = matRad_scanDicomImportFolder(obj)
%
% input
%   patDir:         folder to be scanned
%
% output
%   allfiles:       matlab struct with a list of dicom files including meta
%                   infomation (type, series number etc.)
%   patientList:    list of patients with dicom data in the folder
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

%% print current status of the import script
matRad_cfg.dispInfo('Dose series matched to the different plans are displayed and could be selected.\n');
matRad_cfg.dispInfo('Rechecking of correct matching procedure is recommended.\n');

global warnDlgDICOMtagShown;
warnDlgDICOMtagShown = false;

%% get all files in search directory

% dicom import needs image processing toolbox -> check if available
available = matRad_checkEnvDicomRequirements(matRad_cfg.env);

if ~available
    matRad_cfg.dispError('Image processing toolbox / packages not available!');
end

obj.allfiles = matRad_listAllFiles(obj.patDir);

if ~isempty(obj.allfiles)
    %% check for dicom files and differentiate patients, types, and series
    
    % Arrays of slice locations to find z resolution 
    % if it is not given initially in the files 
    LocationsArray1 = ones(1, numel(obj.allfiles(:,1))); 
    LocationsArray = LocationsArray1*1000;
    ThBool = [];
    numOfFiles = numel(obj.allfiles(:,1));
    h = waitbar(0,'Please wait...','Color',matRad_cfg.gui.backgroundColor,'DefaultTextColor',matRad_cfg.gui.textColor);
    matRad_applyThemeToWaitbar(h);
    % precision value for double to string conversion
    str2numPrc = 10;
    %h.WindowStyle = 'Modal';
    steps = numOfFiles;
    for i = numOfFiles:-1:1
        waitbar((numOfFiles+1-i) / steps)
        try % try to get DicomInfo
            if matRad_cfg.isOctave || verLessThan('matlab','9')
                info = dicominfo(obj.allfiles{i});
            else
                info = dicominfo(obj.allfiles{i},'UseDictionaryVR',true);
            end
        catch
            obj.allfiles(i,:) = [];
            
            % Show progress
            if matRad_cfg.logLevel > 2
                matRad_progress(numOfFiles+1-i, numOfFiles);
            end
            
            continue;
        end
        try
            obj.allfiles{i,2} = info.Modality;
        catch
            obj.allfiles{i,2} = NaN;
        end
        
        obj.allfiles = parseDicomTag(obj.allfiles,info,'PatientID',i,3);
        
        switch obj.allfiles{i,2}
            case 'CT'
               
               obj.allfiles = parseDicomTag(obj.allfiles,info,'SeriesInstanceUID',i,4);
                
            case {'RTPLAN','RTDOSE','RTSTRUCT'}
               
               obj.allfiles = parseDicomTag(obj.allfiles,info,'SOPInstanceUID',i,4);

           otherwise
               
              obj.allfiles = parseDicomTag(obj.allfiles,info,'SeriesInstanceUID',i,4);

        end
        
        obj.allfiles = parseDicomTag(obj.allfiles,info,'SeriesNumber',i,5,@seriesnum2str); %We want to make sure the series number is stored as string
        obj.allfiles = parseDicomTag(obj.allfiles,info,'FamilyName',i,6);
        obj.allfiles = parseDicomTag(obj.allfiles,info,'GivenName',i,7);
        obj.allfiles = parseDicomTag(obj.allfiles,info,'PatientBirthDate',i,8);

        try
            if strcmp(info.Modality,'CT')
                obj.allfiles{i,9} = num2str(info.PixelSpacing(1),str2numPrc);
            else
                obj.allfiles{i,9} = NaN;
            end
        catch
            obj.allfiles{i,9} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')

                obj.allfiles{i,10} = num2str(info.PixelSpacing(2),str2numPrc);
            else
                obj.allfiles{i,10} = NaN;
            end
        catch
            obj.allfiles{i,10} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                %usually the Attribute should be SliceThickness, but it
                %seems like some data uses "SpacingBetweenSlices" instead,
                %but if there is neither this nor that attribute,
                %resolution will be calculated based on SliceLocations
                if isfield(info,'SliceThickness') && info.SliceThickness ~= 0
                    obj.allfiles{i,11} = num2str(info.SliceThickness,str2numPrc);
                    ThBool = 1;
                elseif isfield(info,'SpacingBetweenSlices') && ~isempty(info.SpacingBetweenSlices)
                    obj.allfiles{i,11} = num2str(info.SpacingBetweenSlices,str2numPrc);
                    ThBool = 1;
                else
                    LocationsArray(i) = info.SliceLocation; 
                    
                end
            else
                obj.allfiles{i,11} = NaN;
            end
        catch
            obj.allfiles{i,11} = NaN;
        end
        try

            if strcmp(info.Modality,'RTDOSE')
                dosetext_helper = strcat('Instance','_', num2str(info.InstanceNumber),'_', ...
                    info.DoseSummationType, '_', info.DoseType);
                obj.allfiles{i,12} = dosetext_helper;
            else
                obj.allfiles{i,12} = NaN;
            end
        catch
            obj.allfiles{i,12} = NaN;
        end
        % writing corresponding dose dist.
        try
            if strcmp(obj.allfiles{i,2},'RTPLAN')
                corrDose = [];
                numDose = length(fieldnames(info.ReferencedDoseSequence));
                for j = 1:numDose
                    fieldName = strcat('Item_',num2str(j));
                    corrDose{j} = info.ReferencedDoseSequence.(fieldName).ReferencedSOPInstanceUID;
                end
                obj.allfiles{i,13} = corrDose;
            else
                obj.allfiles{i,13} = {'NaN'};
            end

        catch
            obj.allfiles{i,13} = {'NaN'};
        end
        
        % Show progress
        if matRad_cfg.logLevel > 2
            matRad_progress(numOfFiles+1-i, numOfFiles);
        end
        
    end

    % Filtration, getting and assigning z resolution to all CT files
    FiltredLocArray = unique(LocationsArray);
    FiltredLocArray(end) = [];
    numOfFiles = numel(obj.allfiles(:,1));

    if isempty(ThBool)
        for i = numOfFiles:-1:1
            if strcmp(obj.allfiles{i,2},'CT') 
                obj.allfiles{i,11} = num2str(unique(diff(FiltredLocArray)));
            end
        end
    end    
    
    close(h)
    
    if ~isempty(obj.allfiles)
        obj.patients = unique(obj.allfiles(:,3));
        
        if isempty(obj.patients)
            msgbox('No patient found with DICOM CT _and_ RT structure set in patient directory!', 'Error','error');
        end
    else
        msgbox('No DICOM files found in patient directory!', 'Error','error');
        %h.WindowStyle = 'Modal';
        %error('No DICOM files found in patient directory');
    end
else
    msgbox('Search folder empty!', 'Error','error');
    
end

clear warnDlgDICOMtagShown;

end

function allfiles = parseDicomTag(allfiles,info,tag,row,column,parsefcn)

global warnDlgDICOMtagShown;

defaultPlaceHolder = '001';

if nargin < 6
    parsefcn = @(x) x;
end

try
   if isfield(info,tag)
      if ~isempty(info.(tag))
         allfiles{row,column} = parsefcn(info.(tag));
      else
         allfiles{row,column} = defaultPlaceHolder;
      end
   else
      allfiles{row,column} = defaultPlaceHolder;
   end
catch
   allfiles{row,column} = NaN;
end

if ~warnDlgDICOMtagShown && strcmp(allfiles{row,column},defaultPlaceHolder) && (column == 3 || column == 4)
 
   dlgTitle    = 'Dicom Tag import';
   dlgQuestion = ['matRad_scanDicomImportFolder: Could not parse dicom tag: ' tag '. Using placeholder ' defaultPlaceHolder ' instead. Please check imported data carefully! Do you want to continue?'];
   answer      = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
   
   warnDlgDICOMtagShown = true;
   
   switch answer
      case 'No'
         matRad_cfg.dispError('Inconsistency in DICOM tags')  
   end
end

end

function value = seriesnum2str(value)
    if isnumeric(value)
        value = num2str(value);
    end
end




