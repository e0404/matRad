classdef matRad_DicomImporter < handle
    % matRad_DicomImporter matRad class to handle a dicom import.  
    %
    % Example on how to use the matRad_DicomImport class
    %
    % dcmImpObj = matRad_DicomImporter('pathToFolder');          % create instance of matRad_DicomImporter
    % dcmImpObj.matRad_importDicom(dcmImpObj);                   % run the import
    %
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
       
       % path to DICOM file 
       patDir;

       % lists of all files
       allfiles;
       patients;
       selectedPatient;
       importFiles; % all the names (directories) of files, that will be imported

       % properties with data for import functions
       importCT;
       importRtss;
       importRTDose;
       
       % structures for .mat file
       ct = [];
       cst = [];
       stf = [];
       pln = [];
       resultGUI = [];
       
       ImportGrid;
       
       % bools 
       dicomMetaBool;
       visBool;



    end
    
    methods

        function obj = matRad_DicomImporter(pathToFolder)

            %matRad_DicomImporter Construct an instance of this class
            %   Can be called with the structures. If no argument is given,
            %   all structures will be read from the base workspace

            obj.patDir = pathToFolder;
            matRad_cfg = MatRad_Config.instance();
            
            if matRad_cfg.isOctave
                %Octave needs the DICOM package
                try
                    pkg load dicom;
                catch
                    matRad_cfg.dispError('The DICOM export requires the octave-forge package "dicom"!\n');
                end
            end

            obj.patDir = pathToFolder;
                        
            obj.matRad_scanDicomImportFolder();

            % matRad_DicomImporter imports only one structure, to select
            % patients and structures within a single patient the
            % matRad_importDicomWidget is used

            ctFiles = strcmp(obj.allfiles(:,2),'CT');
            rtssFiles = strcmpi(obj.allfiles(:,2),'rtstruct'); %note we can have multiple RT structure sets, matRad will always import the firstit finds
            rtPlanFiles = strcmpi(obj.allfiles(:,2),'rtplan');
            rtDoseFiles = strcmpi(obj.allfiles(:,2),'rtdose');

            obj.importFiles.ct = obj.allfiles(ctFiles,:);%All CT slice filepaths stored in a cell array like {'CTSlice1.dcm','CTSlice2.dcm'};
            obj.importFiles.rtss = obj.allfiles(rtssFiles,:); 
            obj.importFiles.rtplan = obj.allfiles(rtPlanFiles,:);
            obj.importFiles.rtdose = obj.allfiles(rtDoseFiles,:);
            
            for i = numel(obj.allfiles(:,1)):-1:1
                if strcmp(obj.allfiles(i,2),'CT')
                    obj.importFiles.resx = obj.allfiles{i,9}; 
                    obj.importFiles.resy = obj.allfiles{i,10};
                    obj.importFiles.resz = obj.allfiles{i,11}; %some CT dicoms do not follow  the standard and use SpacingBetweenSlices 
                    break
                end
            end 
             
            obj.importFiles.useImportGrid = false;
            
                   
        end

        matRad_importDicom(obj)

        obj = matRad_importDicomCt(obj)

        obj = matRad_importDicomRTDose(obj)

        obj = matRad_importDicomRTPlan(obj)

        obj = matRad_importDicomRtss(obj)

        obj = matRad_importDicomSteeringPhotons(obj)

        obj = matRad_importDicomSteeringParticles(obj)

        obj = matRad_scanDicomImportFolder(obj)
        
        obj = matRad_calcHU(obj)

        obj = matRad_createCst(obj)
        
        obj = matRad_dummyCst(obj)

        % matRad_saveImport(obj);
        
    end

end

