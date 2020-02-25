classdef exportGUI_alsKlasse < handle
    
    properties
        guiHandle
    end
    
    methods
        function this = exportGUI_alsKlasse
            this.guiHandle = this.createLayout();
            
            %assign guidata like in guide
            handles = guihandles(this.guiHandle);
            guidata(this.guiHandle,handles);
            
            this.initialize();
            
        end
    end
    
    methods
        function guiHandle = createLayout(this)
            h1 = figure(...
                'IntegerHandle','off',...
                'Renderer', 'painters',...
                'MenuBar','none',...
                'NumberTitle','off',...
                'PaperUnits','inches',...
                'Position', [450 170 480 600],...
                'Color',[0.5 0.5 0.5],...
                'Name','Export Patient',...
                'PaperSize',[8.5 11],...
                'Renderer', 'painters', ...
                'PaperType','usletter');
            
            guiHandle = h1;
            
            %EXPORT BUTTON
            h2 = uicontrol(...
                'BackgroundColor',[0.5 0.5 0.5],...
                'Parent',h1,...
                'Units','normalized',... 
                'String','Export',...
                'UserData',[],...
                'Position',[0.4 0.07 0.2 0.07],... 
                'Tag','btn_export',...
                'Callback',@(hObject,eventdata) this.btn_export_Callback(hObject,eventdata,guidata(hObject)));
            
            %CANCEL BUTTON
            h3 = uicontrol(...
                'Parent',h1,...
                'BackgroundColor',[0.5 0.5 0.5],...
                'Units','normalized',...
                'String','Cancel',...
                'Position',[0.7 0.07 0.2 0.07],...
                'Callback',@(hObject,eventdata) this.btn_cancel_Callback(hObject,eventdata,guidata(hObject)),...
                'Tag','btn_cancel');
            
            %CT CHECKBOX
            h4 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',[0.5 0.5 0.5],...
                'String','CT',...
                'Style','checkbox',...
                'Position',[0.035 0.8 0.7 0.04],...
                'Callback',@(hObject,eventdata) this.checkbox_CT_Callback(hObject,eventdata,guidata(hObject)),...
                'Tag','checkbox_CT');
            
            %Compress CHECKBOX
            h14 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',[0.5 0.5 0.5],...
                'HorizontalAlignment','left',...
                'String','Compress',...
                'Style','checkbox',...
                'Value',1,...
                'Position',[0.035 0.15 0.7 0.04],...
                'Tag','checkbox_compress',...
                'Callback',@(hObject,eventdata) this.checkbox_Compress_Callback(hObject,eventdata,guidata(hObject)));
            
           
            % RESULT CUBES CHECKBOX
            h5 = uicontrol(...
                'Parent',h1,...
                'BackgroundColor',[0.5 0.5 0.5],...
                'HorizontalAlignment','left',...
                'Units','normalized',...
                'String','Result Cubes',...
                'Style','checkbox',...
                'Position',[0.035 0.55 0.7 0.04],...
                'Callback',@(hObject,eventdata) this.checkbox_dose_Callback(hObject,eventdata,guidata(hObject)),...
                'Tag','checkbox_dose');
            
            %TEXT
            h6 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                'HorizontalAlignment','left',...
                'String','Select export folder:',...
                'Style','text',...
                'BackgroundColor',[0.5 0.5 0.5],...
                'Position',[0.035 0.9 0.7 0.04 ],...
                'Tag','label_dir_export');
            
            %BROWSER PUSHBUTTON
            h7 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'String','Browse',...
                'BackgroundColor',[0.5 0.5 0.5],...
                'Position',[0.75 0.86 0.15 0.05],... 
                 'Callback',@(hObject,eventdata) this.pushbutton_dir_export_browse_Callback(hObject,eventdata,guidata(hObject)),...
                'Tag','pushbutton_dir_export_browse');
            
            %EDIT TEXTFELD
            h8 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'Style','edit',...
                'Position',[0.035 0.865 0.7 0.04 ],...
                'BackgroundColor',[1 1 1],...
                'Callback',@(hObject,eventdata) this.edit_dir_export_Callback(hObject,eventdata,guidata(hObject)),...
                'Tag','edit_dir_export');
            
            %EXTENSION TEXT
            h9 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'BackgroundColor',[0.5 0.5 0.5],...
                'String','Extension',...
                'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                'Style','text',...
                'Position',[0.035 0.25 0.7 0.04],...
                'Tag','text_extension');
            
            % POP-UP MENU
            h10 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'String',{  '*.nrrd'; '*.vtk'; '*.mha' },...
                'Style','popupmenu',...
                'Value',1,...
                'Position',[0.035 0.21 0.75 0.03],... 
                'BackgroundColor',[1 1 1],...
                'Callback',@(hObject,eventdata) this.popupmenu_extension_Callback(hObject,eventdata,guidata(hObject)),...
                'Tag','popupmenu_extension');
            
            h11 = uitable(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',[1 1 1;0.941176470588235 0.941176470588235 0.941176470588235],...
                'ColumnWidth',{  20 278 },...
                'Visible','off',...    
                'Position',[0.1 0.6 0.73 0.2],...
                'ColumnEditable',false,...
                'ColumnFormat',{  'logical' 'char' },...
                'Tag','uitable_vois');
            
            h12 = uitable(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',[1 1 1;0.941176470588235 0.941176470588235 0.941176470588235],...
                'ColumnName',blanks(0),...
                'ColumnWidth',{  20 278 },...
                'RowName',blanks(0),...
                'Position',[0.1 0.32 0.73 0.2],...
                'Visible','off',...
                'ColumnEditable',false,...
                'ColumnFormat',{  'logical' 'char' },...
                'Tag','uitable_doseCubes',...
                'UserData',[]);
%                 'KeyPressFcn',blanks(0),...
%                 'KeyReleaseFcn',blanks(0)
        end
        
        function this = initialize(this)
            
            handles = guidata(this.guiHandle);
            
            %Fills structure export table
            if evalin('base','exist(''cst'',''var'')') == 1
                cst = evalin( 'base', 'cst' );
                tableData = cell(numel(cst(:,2)),2);
                tableData(:,2) = cst(:,2);
                tableData(:,1) = {true};
            else
                tableData = cell(0);
                set(handles.checkbox_CT,'Enable','off');
            end
            set(handles.uitable_vois,'data',tableData);
            
            %Fills result cubes export table
            if evalin('base','exist(''resultGUI'',''var'')')
                result = evalin( 'base', 'resultGUI' );
                cubeNames = fieldnames(result);
                cubeIx = 1;
                for f = 1:numel(cubeNames)
                    if ndims(result.(cubeNames{f})) < 3
                        continue;
                    end
                    cubes{cubeIx} = cubeNames{f};
                    cubeIx = cubeIx + 1;
                end
                numCubes = cubeIx - 1;
                tableData = cell(numCubes,2);
                tableData(:,2) = cubes;
                tableData(:,1) = {true};
            else
                tableData = cell(0);
                set(handles.checkbox_dose,'Enable','off'); %CHANGED CODE!ALTE VERSION: set(handles.checkbox_dose,'Enable','off');
            end
            set(handles.uitable_doseCubes,'data',tableData);
            
            % Update handles structure
            guidata(this.guiHandle, handles);
        end
    end
    
    methods
        %---------------CALLBACK FOR H2 BUTTON EXPORT
        function this = btn_export_Callback(this,  hObject, event, handles)
            exportDir = get(handles.edit_dir_export,'String');
            
            %Sanity check
            if numel(exportDir) == 0
                errordlg('No Export folder selected!');
                return;
            elseif ~exist(exportDir,'dir')
                errordlg(['Folder ' exportDir ' does not exist!']);
                return;
            else
                %Add file separator if necessary
                if exportDir(end) ~= filesep
                    exportDir = [exportDir filesep];
                end
            end
            
            %Get the file extension
            extensionIndex = get(handles.popupmenu_extension,'Value');
            extensions = get(handles.popupmenu_extension,'String');
            extension = extensions{extensionIndex};
            extension = extension(2:end);
            
            saveCT = get(handles.checkbox_CT,'Value');
            saveResults = get(handles.checkbox_dose,'Value');
            
            if (saveCT)
                voiDir = [exportDir '/vois/'];
                if ~exist(voiDir,'dir')
                    if ~mkdir(voiDir)
                        warndlg('Could not create subfolder for VOI masks. Masks will be stored in base folder.');
                        voiDir = exportDir;
                    end
                end
            end
            
            %If we export results, try to create a subdirectory for VOIs
            if (saveResults)
                resultDir = [exportDir '/results/'];
                if ~exist(resultDir,'dir')
                    if ~mkdir(resultDir)
                        warndlg('Could not create subfolder for resulting dose cubes. Cubes will be stored in base folder.');
                        resultDir = exportDir;
                    end
                end
            end
            
            %prepare metadata
            ct = evalin('base','ct');
            
            metadata.resolution = [ct.resolution.x ct.resolution.y ct.resolution.z];
            metadata.compress = get(handles.checkbox_compress,'Value');
            
            %Check if we have position information
            if isfield(ct,'dicomInfo')
                if isfield(ct.dicomInfo,'ImagePositionPatient')
                    metadata.imageOrigin = ct.dicomInfo.ImagePositionPatient;
                    if ~isrow(metadata.imageOrigin)
                        metadata.imageOrigin = transpose(metadata.imageOrigin);
                    end
                end
            end
            
            %This is only for the waitbar to get the number of cubes you wanna save
            numExportCubes = 0;
            if (saveCT)
                if isfield(ct,'cubeHU')
                    numExportCubes = numExportCubes + 1;
                end
                
                if isfield(ct,'cube')
                    numExportCubes = numExportCubes + 1;
                end
                voiNames = get(handles.uitable_vois,'Data');
                voiIndices = find([voiNames{:,1}] == true);
                numExportCubes = numExportCubes + numel(voiIndices);
                
            else
                numExportCubes = 0;
            end
            
            if saveResults
                cubeNames = get(handles.uitable_doseCubes,'data');
                cubeIndices = find([cubeNames{:,1}] == true);
                numExportCubes = numExportCubes + numel(cubeIndices);
            end
            
            %Give an error if nothing was selected
            if numExportCubes == 0
                errordlg('No data was selected for export!');
                return;
            end
            
            currentCube = 0;
            
            hWaitbar = waitbar(0,'Exporting...','WindowStyle', 'modal');
            cleanUp = onCleanup(@() close(hWaitbar));
            
            %CT and Mask export
            if saveCT
                
                if isfield(ct,'cube')
                    %Export the CT (ED suffix to clarify it is not in HU)
                    currentCube = currentCube + 1;
                    waitbar(currentCube/numExportCubes,hWaitbar,['Exporting CT Intensity values (' num2str(currentCube) '/' num2str(numExportCubes) ') ...']);
                    matRad_writeCube(fullfile(exportDir,['CT_ED' extension]),ct.cube{1},'double',metadata);
                end
                
                if isfield(ct,'cubeHU')
                    currentCube = currentCube + 1;
                    waitbar(currentCube/numExportCubes,hWaitbar,['Exporting CT in HU (' num2str(currentCube) '/' num2str(numExportCubes) ') ...']);
                    matRad_writeCube(fullfile(exportDir,['CT_HU' extension]),ct.cubeHU{1},'double',metadata);
                end
                
                %Export VOI masks
                cst = evalin('base','cst');
                
                for voiIx = voiIndices
                    %Waitbar
                    currentCube = currentCube + 1;
                    waitbar(currentCube/numExportCubes,hWaitbar,['Exporting Segmentation Mask (' num2str(currentCube) '/' num2str(numExportCubes) ') ...']);
                    
                    %Get the index list
                    voiRow = find(strcmp(voiNames{voiIx,2},cst(:,2)));
                    voiIndexList = cst{voiRow,4}{1};
                    %Set up the full mask
                    voiMask = zeros(ct.cubeDim);
                    voiMask(voiIndexList) = 1;
                    %Export...
                    matRad_writeCube(fullfile(voiDir,[voiNames{voiIx,2} extension]),voiMask,'uint8',metadata);
                end
                
            end
            
            %Results Export
            if saveResults
                results = evalin('base','resultGUI');
                cubeNames = get(handles.uitable_doseCubes,'data');
                
                for cubeIx = cubeIndices
                    %Export
                    currentCube = currentCube + 1;
                    waitbar(currentCube/numExportCubes,hWaitbar,['Exporting Results (' num2str(currentCube) '/' num2str(numExportCubes) ') ...']);
                    matRad_writeCube(fullfile(resultDir,[cubeNames{cubeIx,2} extension]),results.(cubeNames{cubeIx,2}),'double',metadata);
                end
            end
            
            %close(data.figure1);
            
        end
        
        %------------CALLBACK FOR H3 BUTTON CANCEL
        function this = btn_cancel_Callback(this, hObject, event, guidata)
            close;
            % close(handles.figure1);
        end
        
        %------------CALLBACK FOR H4 BUTTON CHECKBOX CT
        function this = checkbox_CT_Callback(this, hObject, event, guidata)
            saveCT = get(hObject,'Value');
            %Show the VOI-table only if we want to save a CT
            if (saveCT)
                set(guidata.uitable_vois,'Visible', 'on', 'Enable','on');
            else
                set(guidata.uitable_vois,'Visible', 'off', 'Enable','off');
            end
        end
        
        %-------------CALLBACK FOR H5 BUTTON CHECKBOX DOSE
        function this = checkbox_dose_Callback(this, hObject, event,guidata)
            saveDose = get(hObject,'Value');
            if (saveDose)
                set(guidata.uitable_doseCubes,'Visible', 'on', 'Enable','on');
                
            else
                set(guidata.uitable_doseCubes,'Visible', 'off', 'Enable','off');
            end
        end
        
        %-------------CALLBACK FOR H7 PUSHBUTTON DIR EXPORT BROWSE
        % CAN'T CHANGE THE VARIABLES BECAUSE OF GUIDATA LINE 415
        function this = pushbutton_dir_export_browse_Callback(this, hObject, event,data)
            exportDir = uigetdir('', 'Choose the export directory...');
            if exportDir ~= 0
                exportDir = [exportDir filesep];
                set(data.edit_dir_export,'String',exportDir);
                % Update handles structure
                guidata(hObject, data);
            end
        end
        
        %------------CALLBACK FOR H8 EDIT DIR EXPORT
        % SAME PROBLEM HERE!?
        function this = edit_dir_export_Callback(this, hObject, event,handles)
            exportDir = get(handles.edit_dir_export,'String');
            
            %Add filesperator
            if exportDir(end) ~= filesep
                exportDir = [exportDir filesep];
            end
            
            %Check if the user specified an existing directory
            if ~exist(exportDir,'dir')
                warndlg(['Folder ' exportDir ' does not exist!']);
                exportDir = '';
            end
            
            set(handles.edit_dir_export,'String',exportDir);
            guidata(hObject, handles);
            %guidata(src,data);
        end
        
        %------------CALLBACK FOR H10 POP-UP MENU EXTENSION
        function this = popupmenu_extension_Callback(this, hObject, event,guidata)
            
        end
        
        %------------CALLBACK FOR H14 CHECKBOX COMPRESS
        function this = checkbox_Compress_Callback(this, hObject, event,guidata)

        end
    end
end

