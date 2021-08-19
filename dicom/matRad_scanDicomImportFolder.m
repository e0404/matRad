function [ fileList, patientList ] = matRad_scanDicomImportFolder( patDir )
% matRad function to scan a folder for dicom data
% 
% call
%   [ fileList, patientList ] = matRad_scanDicomImportFolder( patDir )
%
% input
%   patDir:         folder to be scanned
%
% output
%   fileList:       matlab struct with a list of dicom files including meta
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
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
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
if ~license('checkout','image_toolbox')
    matRad_cfg.dispError('image processing toolbox and/or corresponding licence not available');
end

fileList = matRad_listAllFiles(patDir);

if ~isempty(fileList)
    %% check for dicom files and differentiate patients, types, and series
    numOfFiles = numel(fileList(:,1));
    h = waitbar(0,'Please wait...');
    % precision value for double to string conversion
    str2numPrc = 10;
    %h.WindowStyle = 'Modal';
    steps = numOfFiles;
    for i = numOfFiles:-1:1
        waitbar((numOfFiles+1-i) / steps)
        try % try to get DicomInfo
            if verLessThan('matlab','9')
                info = dicominfo(fileList{i});
            else
                info = dicominfo(fileList{i},'UseDictionaryVR',true);
            end
        catch
            fileList(i,:) = [];
            matRad_progress(numOfFiles+1-i, numOfFiles);
            continue;
        end
        try
            fileList{i,2} = info.Modality;
        catch
            fileList{i,2} = NaN;
        end
        
        fileList = parseDicomTag(fileList,info,'PatientID',i,3);
        
        switch fileList{i,2}
            case 'CT'
               
               fileList = parseDicomTag(fileList,info,'SeriesInstanceUID',i,4);
                
            case {'RTPLAN','RTDOSE','RTSTRUCT'}
               
               fileList = parseDicomTag(fileList,info,'SOPInstanceUID',i,4);

           otherwise
               
              fileList = parseDicomTag(fileList,info,'SeriesInstanceUID',i,4);

        end
        
        fileList = parseDicomTag(fileList,info,'SeriesNumber',i,5,@seriesnum2str); %We want to make sure the series number is stored as string
        fileList = parseDicomTag(fileList,info,'FamilyName',i,6);
        fileList = parseDicomTag(fileList,info,'GivenName',i,7);
        fileList = parseDicomTag(fileList,info,'PatientBirthDate',i,8);

        try
            if strcmp(info.Modality,'CT')
                fileList{i,9} = num2str(info.PixelSpacing(1),str2numPrc);
            else
                fileList{i,9} = NaN;
            end
        catch
            fileList{i,9} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                fileList{i,10} = num2str(info.PixelSpacing(2),str2numPrc);
            else
                fileList{i,10} = NaN;
            end
        catch
            fileList{i,10} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                %usually the Attribute should be SliceThickness, but it
                %seems like some data uses "SpacingBetweenSlices" instead.
                if isfield(info,'SliceThickness') && ~isempty(info.SliceThickness)
                    fileList{i,11} = num2str(info.SliceThickness,str2numPrc);
                elseif isfield(info,'SpacingBetweenSlices')
                    fileList{i,11} = num2str(info.SpacingBetweenSlices,str2numPrc);
                else
                    matRad_cfg.dispError('Could not identify spacing between slices since neither ''SliceThickness'' nor ''SpacingBetweenSlices'' are specified');
                end
            else
                fileList{i,11} = NaN;
            end
        catch
            fileList{i,11} = NaN;
        end
        try
            if strcmp(info.Modality,'RTDOSE')
                dosetext_helper = strcat('Instance','_', num2str(info.InstanceNumber),'_', ...
                    info.DoseSummationType, '_', info.DoseType);
                fileList{i,12} = dosetext_helper;
            else
                fileList{i,12} = NaN;
            end
        catch
            fileList{i,12} = NaN;
        end
        % writing corresponding dose dist.
        try
            if strcmp(fileList{i,2},'RTPLAN')
                corrDose = [];
                numDose = length(fieldnames(info.ReferencedDoseSequence));
                for j = 1:numDose
                    fieldName = strcat('Item_',num2str(j));
                    corrDose{j} = info.ReferencedDoseSequence.(fieldName).ReferencedSOPInstanceUID;
                end
                fileList{i,13} = corrDose;
            else
                fileList{i,13} = {'NaN'};
            end

        catch
            fileList{i,13} = {'NaN'};
        end
        matRad_progress(numOfFiles+1-i, numOfFiles);
        
    end
    close(h)
    
    if ~isempty(fileList)
        patientList = unique(fileList(:,3));
        
        if isempty(patientList)
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

function fileList = parseDicomTag(fileList,info,tag,row,column,parsefcn)

global warnDlgDICOMtagShown;

defaultPlaceHolder = '001';

if nargin < 6
    parsefcn = @(x) x;
end

try
   if isfield(info,tag)
      if ~isempty(info.(tag))
         fileList{row,column} = parsefcn(info.(tag));
      else
         fileList{row,column} = defaultPlaceHolder;
      end
   else
      fileList{row,column} = defaultPlaceHolder;
   end
catch
   fileList{row,column} = NaN;
end

if ~warnDlgDICOMtagShown && strcmp(fileList{row,column},defaultPlaceHolder) && (column == 3 || column == 4)
 
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




