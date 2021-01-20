classdef matRad_DicomExporter < handle
    % matRad_DicomExporter matRad class to handle a dicom export.
    %
    %
    % Example on how to use the matRad_DicomExport class
    %
    % dcmExpObj          = matRad_DicomExporter;   % create instance of matRad_DicomExporter
    % dcmExpObj.dicomDir = 'pathAsString';         % set the output path for the Dicom export
    % dcmExp.matRad_exportDicom();                 % run the export
    %
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
       
        %output folder
        dicomDir = ['.' filesep];
        
        % matRad structures to export
        ct = [];
        cst = [];
        stf = [];
        pln = [];
        resultGUI = [];
        
        % Study & Patient Dicom Information
        PatientID
        PatientName
        PatientPosition
        PatientBirthDate = '';
        PatientSex = '';
        StudyID
        StudyDate
        StudyTime
        StudyInstanceUID
        FrameOfReferenceUID
        
        %Operator
        OperatorsName
        
        % ct
        ctFilePrefix = 'ct_slice_';
        ctMeta
        ctSliceMetas
        ctExportStatus
        
        % rtstruct
        rtssFileName = 'RTStruct';
        rtssMeta
        rtssExportStatus
        
        % RTdose
        rtDoseFilePrefix = 'RTDose_';
        rtDoseNames
        rtDoseMetas
        rtDoseExportStatus
        
        % some dictionaries
        externalContourDict = {'EXTERNAL','BODY','PATIENT'}; %Names to identify external contours
        targetPtvDict = {'PTV','MARGIN'};
        targetCtvDict = {'CTV'};
        targetGtvDict = {'GTV','TUMOR'};
    end
    
    methods
        function obj = matRad_DicomExporter(ct,cst,pln,stf,resultGUI)
            %matRad_DicomExporter Construct an instance of this class
            %   Can be called with the structures. If no argument is given,
            %   all structures will be read from the base workspace
            
            
            matRad_cfg = MatRad_Config.instance();  
            
            env = matRad_getEnvironment();
            if strcmp(env,'OCTAVE')
                %Octave needs the DICOM package
                try
                    pkg load dicom;
                catch
                    matRad_cfg.dispError('The DICOM export requires the octave-forge package "dicom"!\n');
                end
            end
            
            if nargin == 0
                try
                    obj.ct = evalin('base','ct');
                catch
                    matRad_cfg.dispInfo('matRad_DicomExporter: No CT found in workspace.\n');
                end
                try
                    obj.cst = evalin('base','cst');
                catch
                    matRad_cfg.dispInfo('matRad_DicomExporter: No cst found in workspace.\n');
                end
                try
                    obj.pln = evalin('base','pln');
                catch
                    matRad_cfg.dispDebug('matRad_DicomExporter: No pln found in workspace.\n');
                end
                try
                    obj.stf = evalin('base','stf');
                catch
                    matRad_cfg.dispDebug('matRad_DicomExporter: No stf found in workspace.\n');
                end
                try
                    obj.resultGUI = evalin('base','resultGUI');
                catch
                end
            else
                if exist('ct','var')
                    obj.ct = ct;
                end
                if exist('cst','var')
                    obj.cst = cst;
                end
                if exist('pln','var')
                    obj.pln = pln;
                end
                if exist('stf','var')
                    obj.stf = stf;
                end
                if exist('resultGUI','var')
                    obj.pln = resultGUI;
                end
            end
            
            obj = obj.generateDefaultDicomData();
        end
        
        function obj = generateDefaultDicomData(obj)
            % generateDefaultDicomData Fill Patient & Study dicom info
                        
            date    = now;
            dateStr = datestr(date,'yyyymmdd');
            timeStr = datestr(date,'HHMMSS');
            
            obj.PatientID = ['matRad_default_patient_' dateStr];
            obj.PatientName = obj.dicomName();
            obj.OperatorsName = obj.dicomName();
            
            obj.StudyID = ['matRad_' dateStr];
            obj.StudyDate = dateStr;
            obj.StudyTime = timeStr;
            obj.StudyInstanceUID = dicomuid;
            
            % coordinates
            obj.FrameOfReferenceUID = dicomuid;
            
            if isfield(obj.ct,'dicomInfo') && isfield(obj.ct.dicomInfo,'PatientPosition')
               obj.PatientPosition = obj.ct.dicomInfo.PatientPosition;
            else
               obj.PatientPosition = 'HFS'; %matRad default coordinate system
            end
        end
        
        obj = matRad_exportDicom(obj)
        
        obj = matRad_exportDicomCt(obj)
        
        obj = matRad_exportDicomRTStruct(obj)
        
        obj = matRad_exportDicomRTPlan(obj)
        
        obj = matRad_exportDicomRTDoses(obj)
        
    end
    
    methods (Static)
        function sarray = addStruct2StructArray(sarray,s,i)
            %addStruct2StructArray Helper function to assign structs
            %   When assigning structs to a struct array, there is usually
            %   an error. This method works around this error by assigning
            %   a struct s to the structure array sarray field by field at
            %   position i.
            %
            %call
            %   sarray = matRad_DicomExporter.matRad_exportDicomRTStruct(sarray,s,i)
            %
            %Input:
            %   sarray: structure array to add struct s to
            %   s:      struct to add to sarray
            %   i:      optional index. if not given, s is appended
            
            if nargin < 3
                i = numel(sarray)+1;
            end
            fnames = fieldnames(s);
            sarray(i).(fnames{1}) = s.(fnames{1});
            if length(fnames) > 1
                for k = 2:length(fnames)
                    sarray(i).(fnames{k}) = s.(fnames{k});
                end
            end
        end
                
        function meta = assignDefaultMetaValue(meta,fn,default,displayBool)
            %assignDefaultMetaValue Helper function to default meta values
            %   When meta information is already given, it may not be
            %   overwritten by a default value. This function automatically
            %   checks for a value and optionally prints to console
            %
            %call
            %   meta = matRad_DicomExporter.assignDefaultMetaValue(meta,fn,default,displayBool)
            %
            %Input:
            %   meta:           structure with dicom meta information
            %   fn:             meta fieldname
            %   default:        default value to assign
            %   displayBool:    optional request for console output.
            %
            %Output:
            %   meta:           meta struct with default value or value
            %                   that was already there before
            
            if nargin < 4
                displayBool = true;
            end
            
            if ~isfield(meta,fn)
                meta.(fn) = default;
                
                if isnumeric(default)
                    default = num2str(default);
                elseif isstruct(default)
                    default = '';
                end
                
                if displayBool
                    fprintf('No value set for ''%s'', using default value ''%s''\n',fn,default);
                end
            end
        end
        
        function name = dicomName(family,given,middle,prefix,suffix)
            if nargin < 5
                suffix = '';
            end
            if nargin < 4
                prefix = '';
            end
            if nargin < 3
                middle = '';
            end
            if nargin < 2
                given = '';
            end
            if nargin < 1
                family = '';
            end
            
            name.FamilyName = family;
            name.GivenName = given;
            name.MiddleName = middle;
            name.NamePrefix = prefix;
            name.NameSuffix = suffix;
            
            env = matRad_getEnvironment();
            if strcmp(env,'OCTAVE')
                name = strjoin(struct2cell(name));
            end
            
        end
    end
    
    
end

