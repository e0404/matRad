classdef matRad_importDicomWidget < matRad_Widget
    
    properties
    end
    
    methods
        
        function this = matRad_importDicomWidget(handleParent)
            %MATRAD_IMPORTDICOMWIDGET Construct an instance of this class
            if nargin < 1
                handleParent = figure(...
                    'IntegerHandle','off',...
                    'Renderer', 'painters',...
                    'Units','characters',...
                    'MenuBar','none',...
                    'NumberTitle','off',...
                    'PaperUnits','inches',...
                    'Position', [450 170 480 600],...
                    'Color',[0.5 0.5 0.5],...
                    'Name','Import Dicom',...
                    'PaperSize',[8.5 11],...
                    'Position',[135.8 21.6153846153846 219.4 39.6923076923077],...
                    'Renderer', 'painters', ...
                    'PaperType','usletter');
            end
            this = this@matRad_Widget(handleParent);
        end
        
        % INITIALIZE FUNCTION
        function this = initialize(this)
            handles = this.handles;
            handles.output = handles;
            
            axes(handles.axesMatRadLogo)
            [path,name,ext] = fileparts(mfilename('fullpath'));
            
            [im, ~, alpha] = imread([path filesep '..' filesep '..' filesep 'gfx/matrad_logo.png']);
            q = image(im);
            axis equal off
            set(q, 'AlphaData', alpha);
            % show dkfz logo
            axes(handles.axesDKFZLogo)
            [im, ~, alpha] = imread([path filesep '..' filesep  '..' filesep 'gfx/DKFZ_Logo.png']);
            p = image(im);
            axis equal off
            set(p, 'AlphaData', alpha);
            
            % Update handles structure
            % guidata(hObject, handles);
            this.handles = handles;
        end
        
        % OUTPUT FUNCTION
        function varargout = matRad_importDicomGUI_OutputFcn(this, hObject, eventdata)
            handles = this.handles;
            varargout{1} = handles.output;
            this.handles = handles;
        end
        
    end
    
    methods (Access = protected)
        % CREATE FUNCTION
        this = createLayout(this,handleParent);
    end
    
    % CALLBACKS
    methods
        
        % H15 BROWSER BUTTON CALLBACK
        function patDir = browse_button_Callback(this, hObject, eventdata)
            handles = this.handles;
            patDir = uigetdir('', 'Choose the input directory...');
            if patDir ~= 0
                patDir = [patDir filesep];
                %handles.dir_path_field.String = patDir;
                set(handles.dir_path_field,'String',patDir);
                % Update handles structure
                % guidata(hObject, handles);
                this.handles = handles;
                this.scan(hObject, eventdata)
            end
        end
        
        % H17 PATIENT LISTBOX CALLBACK
        function this = patient_listbox_Callback(this, hObject, eventdata)
            handles = this.handles;
            if ~isempty(get(hObject,'String'))
                % enable Import button
                set(handles.import_button,'Enable','on');
                
                % handles.filelist:
                %   1. Filepath
                %   2. Modality
                %   3. PatientID
                %   4. SeriesUID
                %   5. SeriesNumber
                %   9. res_x
                %   10. res_y
                %   11. res_z
                %   12. detailed dose description - currently not in use for GUI user
                patient_listbox = get(handles.patient_listbox,'String');
                selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
                % this gets a list of rtss series for this patient
                set(handles.rtseries_listbox,'Value',1); % set dummy value to one
                set(handles.rtseries_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),4));
                % this gets a list of rt plan series for this patient
                set(handles.rtplan_listbox,'Value',[]); % set dummy value to none
                set(handles.rtplan_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTPLAN') & strcmp(handles.fileList(:,3), selected_patient),4));
                % this gets a list of dose series for this patient
                set(handles.doseseries_listbox,'Value',[]); % set dummy value to none
                set(handles.doseseries_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTDOSE') & strcmp(handles.fileList(:,3), selected_patient),4));
                % selectedDose
                
                if get(handles.SeriesUID_radiobutton,'Value') == 1
                    % this gets a list of ct series for this patient
                    set(handles.ctseries_listbox,'Value',1); % set dummy value to one
                    set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),4)));
                    
                    selectedDoseSeriesString = get(handles.doseseries_listbox,'String');
                    % this gets a resolution for this patient
                    selectedCtSeriesString = get(handles.ctseries_listbox,'String');
                    if ~isempty(selectedCtSeriesString)
                        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
                        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
                        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
                    else
                        res_x = NaN; res_y = NaN; res_z = NaN;
                    end
                else
                    set(handles.ctseries_listbox,'Value',1); % set dummy value to one
                    set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),5)));
                    selectedCtSeriesString = get(handles.ctseries_listbox,'String');
                    if ~isempty(selectedCtSeriesString)
                        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
                        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
                        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
                    else
                        res_x = NaN; res_y = NaN; res_z = NaN;
                    end
                end
                set(handles.resx_edit,'String',res_x);
                set(handles.resy_edit,'String',res_y);
                if numel(res_z) > 1
                    set(handles.resz_edit,'String','not equi');
                else
                    set(handles.resz_edit,'String',res_z);
                end
                % Update handles structure
                % guidata(hObject, handles);
                this.handles = handles;
            end
        end      
        
        % H22 IMPORT BUTTON CALLBACK
        function this = import_button_Callback(this, hObject, eventdata)
            handles = this.handles;
            patient_listbox = get(handles.patient_listbox,'String');
            ctseries_listbox = get(handles.ctseries_listbox,'String');
            rtplan_listbox = get(handles.rtplan_listbox,'String');
            doseseries_listbox = get(handles.rtplan_listbox,'String');
            if ~isempty(patient_listbox)
                selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
            end
            if ~isempty(ctseries_listbox)
                selected_ctseries = ctseries_listbox(get(handles.ctseries_listbox,'Value'));
            end
            if ~isempty(rtplan_listbox)
                selected_rtplan = rtplan_listbox(get(handles.rtplan_listbox,'Value'));
            end
            
            if get(handles.SeriesUID_radiobutton,'Value') == 1
                files.ct = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
                    strcmp(handles.fileList(:,4), selected_ctseries),:);
                
                %files.rtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
                %    strcmp(handles.fileList(:,4), selected_rtseries),:);
            else
                files.ct = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
                    strcmp(cellfun(@num2str,handles.fileList(:,5),'UniformOutput',false), selected_ctseries) & strcmp(handles.fileList(:,2),'CT'),:);
                
                %files.rtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
                %    strcmp(handles.fileList(:,5), selected_rtseries),:);
            end
            
            allRtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,2),'RTSTRUCT'),:);
            if ~isempty(allRtss)
                files.rtss = allRtss(get(handles.rtseries_listbox,'Value'),:);
            else
                files.rtss = [];
            end
            
            files.resx = str2double(get(handles.resx_edit,'String'));
            files.resy = str2double(get(handles.resy_edit,'String'));
            % check if valid assignment is made when z slices are not equi-distant
            if strcmp(get(handles.resz_edit,'String'),'not equi')
                msgbox('Ct data not sliced equi-distantly in z direction! Chose uniform resolution.', 'Error','error');
                return;
            else
                files.resz = str2double(get(handles.resz_edit,'String'));
            end
            % selected RT Plan
            rtplan = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,2),'RTPLAN'),:);
            if ~isempty(rtplan) && ~isempty(get(handles.rtplan_listbox,'Value'))
                files.rtplan = rtplan(get(handles.rtplan_listbox,'Value'),:);
            end
            
            % selected RT Dose
            rtdose = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,2),'RTDOSE'),:);
            if ~isempty(rtdose) && ~isempty(get(handles.doseseries_listbox,'Value'))
                selectedRtDose   = get(handles.doseseries_listbox,'String');
                selectedRtDoseIx = NaN*ones(1,numel(selectedRtDose));
                for i = 1:numel(selectedRtDose)
                    selectedRtDoseIx(i) = find(strcmp(rtdose(:,4),selectedRtDose{i}));
                end
                files.rtdose = rtdose(selectedRtDoseIx,:);
            end
            
            % check if we should use the dose grid resolution
            files.useDoseGrid = get(handles.checkbox3,'Value');
            
            % dicomMetaBool: store complete DICOM information and patientName or not
            dicomMetaBool = logical(get(handles.checkPatientName,'Value'));
            matRad_importDicom(files, dicomMetaBool);
            
            this.handles = handles;
        end
        
        % H23 CANCEL BUTTON CALLBACK
        function this = cancel_button_Callback(this, hObject, eventdata)
            
            % handles = this.handles;
            % close(handles.figure1);
            %Maybe needs adaptation for use in figures
            close;
            % this.handles = handles;
        end
        
        % H?? RESCAN BUTTON CALLBACK; NICHT IN CREATELAYOUT VORHANDEN
        function this = rescan_button_Callback(this, hObject, eventdata)
        end
        
        % H24 SERIES UID RADIOBUTTON CALLBACK
        function this = SeriesUID_radiobutton_Callback(this, hObject, eventdata)
            handles = this.handles;
            if get(hObject,'Value') == 1
                set(handles.SeriesNumber_radiobutton,'Value',0);
            else
                set(hObject,'Value',1);
                set(handles.SeriesNumber_radiobutton,'Value',0);
            end
            if isfield(handles, 'fileList')
                patient_listbox = get(handles.patient_listbox,'String');
                selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
                set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),4)));
                set(handles.rtseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),4)));
                set(handles.doseseries_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTDOSE') & strcmp(handles.fileList(:,3), selected_patient),4));
                set(handles.rtplan_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTPLAN') & strcmp(handles.fileList(:,3), selected_patient),4)));
            else
                fprintf('No patient loaded, so just switching default display option to SeriesUID. \n');
            end
            %guidata(hObject, handles);
            this.handles = handles;
        end
        
        % H25 SERIESNUMBER RADIO BUTTON CALLBACK
        function this = SeriesNumber_radiobutton_Callback(this, hObject, eventdata)
            handles = this.handles;
            
            if get(hObject,'Value') == 1
                set(handles.SeriesUID_radiobutton,'Value',0);
            else
                set(hObject,'Value',1);
                set(handles.SeriesUID_radiobutton,'Value',0);
            end
            if isfield(handles, 'fileList')
                patient_listbox = get(handles.patient_listbox,'String');
                selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
               % set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),5)));
                set(handles.ctseries_listbox,'String',unique(cell2mat(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),5))));

            else
                fprintf('No patient loaded, so just switching default display option to SeriesNumber. \n');
            end
            % guidata(hObject, handles);
            this.handles = handles;
        end
        
         
        % H34 DOSESERIES LISTBOX CALLBACK
        function this = doseseries_listbox_Callback(this, hObject, eventdata)
            handles = this.handles;
            
            if ~isempty(get(hObject,'Value'))
                set(handles.checkbox3,'Enable','on');
            else
                set(handles.checkbox3,'Value',0);
                set(handles.checkbox3,'Enable','off');
                % retrieve and display resolution for DICOM ct cube
                patient_listbox = get(handles.patient_listbox,'String');
                selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
                selectedCtSeriesString = get(handles.ctseries_listbox,'String');
                if get(handles.SeriesUID_radiobutton,'Value') == 1
                    if ~isempty(selectedCtSeriesString)
                        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
                        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
                        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
                    else
                        res_x = NaN; res_y = NaN; res_z = NaN;
                    end
                else
                    if ~isempty(selectedCtSeriesString)
                        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
                        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
                        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
                    else
                        res_x = NaN; res_y = NaN; res_z = NaN;
                    end
                end
                set(handles.resx_edit,'String',res_x);
                set(handles.resy_edit,'String',res_y);
                if numel(res_z) > 1
                    set(handles.resz_edit,'String','not equi');
                else
                    set(handles.resz_edit,'String',res_z);
                end
                
            end
            
            this.handles = handles;
        end
        
        % H35 RTPLAN LISTBOX CALLBACK
        function this = rtplan_listbox_Callback(this, hObject, eventdata)
            handles = this.handles;
            
            contents = cellstr(get(hObject,'String'));
            if ~isempty(get(hObject,'Value')) && numel(get(hObject,'Value')) == 1
                
                selectedPlan = contents{get(hObject,'Value')};
                % point at plan in listbox
                selectedPlanLoc = strcmp(handles.fileList(:,4),selectedPlan);
                
                % show only the doses corresponding to the plan
                corrDoses = [handles.fileList{selectedPlanLoc,13}];
                numOfDoses = numel(corrDoses);
                corrDosesLoc = zeros(size(handles.fileList(:,1),1),1);
                for j = 1:numOfDoses
                    if ~isnan(corrDoses{j})
                        corrDosesLoc = corrDosesLoc | strcmp(handles.fileList(:,4),corrDoses{j});
                    end
                end
                
                if sum(corrDosesLoc) == 0
                    warndlg('no rt dose file directly associated to plan file. showing all rt dose files.');
                    corrDosesLoc = strcmp(handles.fileList(:,2),'RTDOSE');
                end
                
                set(handles.doseseries_listbox,'Value',[]); % set dummy value to one
                set(handles.doseseries_listbox,'String',handles.fileList(corrDosesLoc,4));
                
                % disable checkbox for use dose grid is currently checked
                if get(handles.checkbox3,'Value') == 1
                    set(handles.checkbox3,'Value',0);
                    checkbox3_Callback(handles.checkbox3,[], handles);
                end
                set(handles.checkbox3,'Enable','off');
                
            elseif numel(get(hObject,'Value')) >=2
                warning('More than one RTPLAN selected. Unsetting selection ...');
                patient_listbox_Callback(this, hObject, eventdata);
            else
                patient_listbox_Callback(this, hObject, eventdata);
            end
            
            this.handles = handles;
            
        end
        
        % H36 DIR PATH FIELD CALLBACK
        function this = dir_path_field_Callback(this, hObject, eventdata)
            handles = this.handles;
            patDir = get(handles.dir_path_field,'String');
            if patDir(end) ~= filesep;
                patDir = [patDir filesep];
                set(handles.dir_path_field,'String',patDir);
                % guidata(hObject, handles);
                this.handles = handles;
            end
            scan(hObject, eventdata);
        end
        
        % H37 CHECK PATIENTNAME CALLBACK
        function this = checkPatientName_Callback(this, hObject, eventdata)
            handles = this.handles;
            % hObject    handle to checkPatientName (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            %A = get(hObject,'Value');
            
            % Hint: get(hObject,'Value') returns toggle state of checkPatientName
            %guidata(hObject, handles);
            this.handles = handles;
            
        end
        
        % H?? CHECKBOX§ CALLBACK
        function this = checkbox3_Callback(this, hObject, eventdata)
            % hObject    handle to checkbox3 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hint: get(hObject,'Value') returns toggle state of checkbox3
            handles = this.handles;
            
            if get(hObject,'Value')
                set(handles.resx_edit,'Enable', 'off');
                set(handles.resy_edit,'Enable', 'off');
                set(handles.resz_edit,'Enable', 'off');
                % retrieve and display resolution for DICOM dose cube
                doseFilesInList = get(handles.doseseries_listbox,'String');
                selectedDoseFiles = get(handles.doseseries_listbox,'Value');
                if isempty(selectedDoseFiles)
                    set(hObject,'Value',0)
                    errordlg('no dose file selected');
                    return;
                end
                for i = 1:numel(selectedDoseFiles)
                    selectedDoseFile = doseFilesInList{selectedDoseFiles(i)};
                    if verLessThan('matlab','9')
                        dicomDoseInfo = dicominfo(handles.fileList{find(strcmp(handles.fileList(:,4),selectedDoseFile)),1});
                    else
                        dicomDoseInfo = dicominfo(handles.fileList{find(strcmp(handles.fileList(:,4),selectedDoseFile)),1},'UseDictionaryVR',true);
                    end
                    res_x{i} = dicomDoseInfo.PixelSpacing(1);
                    res_y{i} = dicomDoseInfo.PixelSpacing(2);
                    res_z{i} = dicomDoseInfo.SliceThickness;
                end
                
                if numel(unique(cell2mat(res_x)))*numel(unique(cell2mat(res_y)))*numel(unique(cell2mat(res_z))) ~= 1
                    set(handles.checkbox3,'Value',0);
                    warndlg('Different resolutions in dose file(s)');
                    set(handles.resx_edit,'Enable', 'on');
                    set(handles.resy_edit,'Enable', 'on');
                    set(handles.resz_edit,'Enable', 'on');
                else
                    set(handles.resx_edit,'String',num2str(res_x{1}));
                    set(handles.resy_edit,'String',num2str(res_y{1}));
                    set(handles.resz_edit,'String',num2str(res_z{1}));
                end
                
            else
                set(handles.resx_edit,'Enable', 'on');
                set(handles.resy_edit,'Enable', 'on');
                set(handles.resz_edit,'Enable', 'on');
                % retrieve and display resolution for DICOM ct cube
                patient_listbox = get(handles.patient_listbox,'String');
                selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
                selectedCtSeriesString = get(handles.ctseries_listbox,'String');
                if get(handles.SeriesUID_radiobutton,'Value') == 1
                    if ~isempty(selectedCtSeriesString)
                        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
                        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
                        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
                    else
                        res_x = NaN; res_y = NaN; res_z = NaN;
                    end
                else
                    if ~isempty(selectedCtSeriesString)
                        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
                        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
                        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
                    else
                        res_x = NaN; res_y = NaN; res_z = NaN;
                    end
                end
                set(handles.resx_edit,'String',res_x);
                set(handles.resy_edit,'String',res_y);
                if numel(res_z) > 1
                    set(handles.resz_edit,'String','not equi');
                else
                    set(handles.resz_edit,'String',res_z);
                end
            end
            
            this.handles = handles;
        end
        
    end
         
        methods (Access = private)
            % SCAN FUNKTION
            function this = scan(this, hObject, eventdata)
                handles = this.handles;
                [fileList, patient_listbox] = matRad_scanDicomImportFolder(get(handles.dir_path_field,'String'));
                if iscell(patient_listbox)
                    handles.fileList =  fileList;
                    %handles.patient_listbox.String = patient_listbox;
                    set(handles.patient_listbox,'String',patient_listbox,'Value',1);
                    % guidata(hObject, handles);
                    this.handles = handles;
                end
            end           
        end
        
    end
    
