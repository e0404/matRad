classdef matRad_DicomExporter < handle
    %matRad_DicomExporter matRad class to handle a dicom export
    
    properties
        %output folder
        dicomDir = ['.' filesep];
        
        %MatRad Structures to export
        ct = [];
        cst = [];
        stf = [];
        pln = [];
        resultGUI = [];
        
        %Study & Patient Dicom Information
        PatientID
        PatientName
        PatientPosition
        StudyID
        StudyDate
        StudyTime
        StudyInstanceUID
        FrameOfReferenceUID
        
        %ct
        ctFilePrefix = 'ct_slice_';
        ctMeta
        ctSliceMetas
        ctExportStatus
        
        %rtstruct
        rtssFileName = 'RTStruct';
        rtssMeta
        rtssExportStatus
        
        %RTdose
        rtDoseFilePrefix = 'RTDose_';
        rtDoseNames
        rtDoseMetas
        rtDoseExportStatus
    end
    
    methods
        function obj = matRad_DicomExporter(ct,cst,pln,stf,resultGUI)
            %matRad_DicomExporter Construct an instance of this class
            %   Can be called with the structures. If no argument is given,
            %   all structures will be read from the base workspace
            if nargin == 0
                try
                    obj.ct = evalin('base','ct');
                end
                try
                    obj.cst = evalin('base','cst');
                end
                try
                    obj.pln = evalin('base','pln');
                end
                try
                    obj.stf = evalin('base','stf');
                end
                try
                    obj.resultGUI = evalin('base','resultGUI');
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
            %generateDefaultDicomData Fill Patient & Study dicom info
            
            date = datetime;
            dateStr = datestr(date,'yyyymmdd');
            timeStr = datestr(date,'HHMMSS');
            
            obj.PatientID = ['matRad_default_patient_' dateStr];
            obj.PatientName = struct('FamilyName','','GivenName','','MiddleName','','NamePrefix','','NameSuffix','');
            obj.StudyID = ['matRad_' dateStr];
            obj.StudyDate = dateStr;
            obj.StudyTime = timeStr;
            obj.StudyInstanceUID = dicomuid;
            
            %Coordinates
            obj.FrameOfReferenceUID = dicomuid;
            obj.PatientPosition = 'HFS'; %matRad default coordinate system
            
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
        
        function dcmName = createEmptyPersonName()
            %createEmptyPersonName Creates a dicom name struct
            dcmName = struct('FamilyName','','GivenName','','MiddleName','','NamePrefix','','NameSuffix','');
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
    end
    
    
end

