classdef matRad_exportWidget < matRad_Widget

    % matRad_exportWidget class to generate GUI widget to export plan as
    % dicom, nrrd etc.
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
    end
    
    methods
        function this = matRad_exportWidget(handleParent)
            matRad_cfg = MatRad_Config.instance();
            if nargin < 1
                handleParent = figure(...
                    'IntegerHandle','off',...
                    'MenuBar','none',...
                    'NumberTitle','off',...
                    'PaperUnits','inches',...
                    'Position', [450 170 440 500],...
                    'Color',matRad_cfg.gui.backgroundColor,...
                    'Name','Export Patient');
            end
            this = this@matRad_Widget(handleParent);
        end
        
        function this = update(this,evt)
            
            doUpdate = true;
            if nargin == 2
                doUpdate = this.checkUpdateNecessary({'resultGUI','ct','cst'},evt);
            end
            
            if ~doUpdate
                return;
            end
            
            handles = this.handles;
            % handles = guidata(this.widgetHandle);
            
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
            this.handles = handles;
        end
    end
    
    
    methods (Access = protected)
        
        function this = createLayout(this)
            
            h1 = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
            
            %EXPORT BUTTON
            h2 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'String','Export',...
                'Tooltip', 'Export selected quantites to selected folder',...
                'UserData',[],...
                'Position',[0.75 0.025 0.2 0.05],...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Tag','btn_export',...
                'Callback',@this.btn_export_Callback);            
            
            %CT CHECKBOX
            h4 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'String','CT',...
                'Tooltip', 'Export CT of selected structures',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Style','checkbox',...
                'Position',[0.035 0.8 0.7 0.05],...
                'Callback',@this.checkbox_CT_Callback,...
                'Tag','checkbox_CT');
            
            
            %Compress CHECKBOX
            h14 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Compress',...
                'Tooltip', 'Export compressed data',...
                'Style','checkbox',...
                'Value',1,...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Position',[0.035 0.13 0.7 0.05],...
                'Tag','checkbox_compress',...
                'Callback',@this.checkbox_Compress_Callback);
            
            
            % RESULT CUBES CHECKBOX
            h5 = uicontrol(...
                'Parent',h1,...
                'HorizontalAlignment','left',...
                'Units','normalized',...
                'String','Result Cubes',...
                'Tooltip', 'Export selected result cubes',...
                'Style','checkbox',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Position',[0.035 0.45 0.7 0.04],...
                'Callback',@this.checkbox_dose_Callback,...
                'Tag','checkbox_dose');
            
            %TEXT
            h6 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Select export folder:',...
                'Tooltip', 'Select the folder you want to export to',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Style','text',...
                'Position',[0.035 0.925 0.915 0.05 ],...
                'Tag','label_dir_export');
            
            %BROWSER PUSHBUTTON
            h7 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'String','Browse',...
                'Tooltip', 'Choose the export directory',...
                'Position',[0.75 0.875 0.2 0.05],...
                'Callback',@this.pushbutton_dir_export_browse_Callback,...
                'Tag','pushbutton_dir_export_browse');
            
            %EDIT TEXTFELD
            h8 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'Style','edit',...
                'Position',[0.035 0.875 0.7 0.05 ],...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Callback',@this.edit_dir_export_Callback,...
                'Tooltip', 'Export path',...
                'Tag','edit_dir_export');
            
            %EXTENSION TEXT
            h9 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Extension',...
                'Tooltip', 'Select file format',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Style','text',...
                'Position',[0.035 0.225 0.7 0.05],...
                'Tag','text_extension');
            
            % DROPDOWN MENU
            h10 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'String',{  '*.nrrd'; '*.vtk'; '*.mha' },...
                'Tooltip', 'File format',...
                'Style','popupmenu',...
                'Value',1,...
                'Position',[0.035 0.21 0.915 0.03],...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Callback',@this.popupmenu_extension_Callback,...
                'Tag','popupmenu_extension');
            
            % CT-TABLE
            h11 = uitable(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'ColumnWidth',{20 400},...
                'Position',[0.035 0.5 0.915 0.3],...
                'ColumnEditable',[true false],...
                'ColumnFormat',{  'logical' 'char' },...
                'ColumnName',{},...
                'RowName',{},...
                'Enable','off',...
                'Tag','uitable_vois');
            
            set(h11,'Units','pixels');
            pos = get(h11,'Position');
            set(h11,'ColumnWidth',{20 pos(3)-39});
            
            
            %RESLT CUBES-TABLES
            h12 = uitable(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'ColumnName',{},...
                'ColumnWidth',{20 400},...
                'RowName',{},...
                'Position',[0.035 0.3 0.915 0.15],...
                'ColumnEditable',[true false],...
                'ColumnFormat',{  'logical' 'char' },...
                'Tag','uitable_doseCubes',...
                'Enable','off',...
                'UserData',[]);
            
            set(h12,'Units','pixels');
            pos = get(h12,'Position');
            set(h12,'ColumnWidth',{20 pos(3)-39});
            
            this.createHandles();
        end
        
    end
    
    
    %CALLBACK METHODAS
    methods
        %---------------CALLBACK FOR H2 BUTTON EXPORT
        function this = btn_export_Callback(this,  hObject, event)
            handles = this.handles;
            
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
                this.handles = handles;
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
            
            try
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
            catch ME
                warning(ME.identifier,'couldn''t export! Reason: %s\n',ME.message)
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
        function this = checkbox_CT_Callback(this, hObject, event)
            handles = this.handles;
            saveCT = get(hObject,'Value');
            %Show the VOI-table only if we want to save a CT
            if (saveCT)
                set(handles.uitable_vois, 'Enable','on');
            else
                set(handles.uitable_vois,'Enable','off');
            end
            this.handles = handles;
        end
        
        %-------------CALLBACK FOR H5 BUTTON CHECKBOX DOSE
        function this = checkbox_dose_Callback(this, hObject, event)
            handles = this.handles;
            saveDose = get(hObject,'Value');
            if (saveDose)
                set(handles.uitable_doseCubes, 'Enable','on');
                %ORIGINAL: set(guidata.uitable_doseCubes,'Visible', 'on', 'Enable','on');
                
            else
                set(handles.uitable_doseCubes,'Enable','off');
                %ORIGINAL: set(guidata.uitable_doseCubes,'Visible', 'off', 'Enable','off');
                
            end
            this.handles = handles;
        end
        
        %-------------CALLBACK FOR H7 PUSHBUTTON DIR EXPORT BROWSE
        function this = pushbutton_dir_export_browse_Callback(this, hObject, event)
            handles = this.handles;
            
            exportDir = uigetdir('', 'Choose the export directory...');
            if exportDir ~= 0
                exportDir = [exportDir filesep];
                set(handles.edit_dir_export,'String',exportDir);
                % Update handles structure
                this.handles = handles;
            end
        end
        
        %------------CALLBACK FOR H8 EDIT DIR EXPORT
        function this = edit_dir_export_Callback(this, hObject, event)
            handles = this.handles;
            
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
            this.handles = handles;
        end
        
        %------------CALLBACK FOR H10 DROPDOWN MENU EXTENSION
        function this = popupmenu_extension_Callback(this, hObject, event)
            
        end
        
        %------------CALLBACK FOR H14 CHECKBOX COMPRESS
        function this = checkbox_Compress_Callback(this, hObject, event)
            
        end
    end
end



