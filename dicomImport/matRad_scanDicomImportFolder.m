function [ fileList, patientList ] = matRad_scanDicomImportFolder( patDir )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% print current status of the import script
fprintf('Dose series matched to the different plans are displayed and could be selected.\n');
fprintf('Rechecking of correct matching procedure is recommended.\n');

%% get all files in search directory

% dicom import needs image processing toolbox -> check if available
if ~license('checkout','image_toolbox')
    error('image processing toolbox and/or corresponding licence not available');
end

fileList = matRad_listAllFiles(patDir);

if ~isempty(fileList)
    %% check for dicom files and differentiate patients, types, and series
    numOfFiles = numel(fileList(:,1));
    h = waitbar(0,'Please wait...');
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
        try
            fileList{i,3} = info.PatientID;
        catch
            fileList{i,3} = NaN;
        end
        switch fileList{i,2}
            case 'CT'
                try
                    fileList{i,4} = info.SeriesInstanceUID;
                catch
                    fileList{i,4} = NaN;
                end
            case 'RTPLAN'
                try
                    fileList{i,4} = info.SOPInstanceUID;
                catch
                    fileList{i,4} = NaN;
                end
            case 'RTDOSE'
                try
                    fileList{i,4} = info.SOPInstanceUID;
                catch
                    fileList{i,4} = NaN;
                end
            case 'RTSTRUCT'
                try
                    fileList{i,4} = info.SOPInstanceUID;
                catch
                    fileList{i,4} = NaN;
               end
            otherwise
                try
                    fileList{i,4} = info.SeriesInstanceUID;
                catch
                    fileList{i,4} = NaN;
                end
        end
        try
            fileList{i,5} = num2str(info.SeriesNumber);
        catch
            fileList{i,5} = NaN;
        end
        try
            fileList{i,6} = info.PatientName.FamilyName;
        catch
            fileList{i,6} = NaN;
        end
        try
            fileList{i,7} = info.PatientName.GivenName;
        catch
            fileList{i,7} = NaN;
        end
        try
            fileList{i,8} = info.PatientBirthDate;
        catch
            fileList{i,8} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                fileList{i,9} = num2str(info.PixelSpacing(1));
            else
                fileList{i,9} = NaN;
            end
        catch
            fileList{i,9} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                fileList{i,10} = num2str(info.PixelSpacing(2));
            else
                fileList{i,10} = NaN;
            end
        catch
            fileList{i,10} = NaN;
        end
        try
            if strcmp(info.Modality,'CT')
                fileList{i,11} = num2str(info.SliceThickness);
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


end

