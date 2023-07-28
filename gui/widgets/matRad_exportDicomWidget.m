classdef matRad_exportDicomWidget < matRad_Widget
    % matRad_exportDicomWidget class to generate GUI widget to export plan
    % to dicom files
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
        variables = {'ct','cst','resultGUI'}; %variables to export
    end
    
    methods
        function this = matRad_exportDicomWidget(handleParent)
            matRad_cfg = MatRad_Config.instance();
            if nargin < 1
                handleParent = figure(...
                    'IntegerHandle','off',...
                    'MenuBar','none',...
                    'NumberTitle','off',...
                    'PaperUnits','inches',...
                    'Position', [450 170 480 350],...
                    'Color',matRad_cfg.gui.backgroundColor,...
                    'Name','Export Dicom');
            end
            this = this@matRad_Widget(handleParent);
            
            update(this);
        end
        function this = initialize(this)
             
        end
        
        
        function this = update(this,evt)
            doUpdate = true;
            if nargin == 2
                %At pln changes and at cst/cst (for Isocenter and new settings) 
                %we need to update
                doUpdate = this.checkUpdateNecessary(this.variables,evt);
            end
            
            if doUpdate                
                handles = this.handles;
                
                % load table with available variables that can be exported
                vExists=false(1,numel(this.variables));
                for i= 1:numel(this.variables)
                    var= char(this.variables(i));
                    vExists(i) = evalin('base',['exist(''' var ''',''var'')']);
                end
                
                if find(vExists,1) % not empty
                    tableData(:,2)= this.variables(vExists);
                    tableData(:,1) ={true};
                    set(handles.uitable_variables,'ColumnEditable',[true,false]);
                else
                    tableData(1,2) = {'No variables to export'};
                    set(handles.btn_export,'Enable','off');
                    set(handles.uitable_variables,'ColumnEditable',[false,false]);
                end
                set(handles.uitable_variables,'data',tableData);
                
                this.handles = handles;
            end
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
                'UserData',[],...
                'Position',[0.75 0.1 0.2 0.1],...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'Tooltip', 'Export selected variables to selected folder',...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Tag','btn_export',...
                'Callback',@this.btn_export_Callback);

            %TEXT
            h6 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                'HorizontalAlignment','left',...
                'String','Select export folder:',...
                'Tooltip', 'Select the folder you want to export to',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Style','text',...
                'Position',[0.035 0.8 0.7 0.1 ],...
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
                'Position',[0.75 0.75 0.2 0.1],...
                'Callback',@this.pushbutton_dir_export_browse_Callback,...
                'Tag','pushbutton_dir_export_browse');
            
            %EDIT TEXTFELD
            h8 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'Style','edit',...
                'Position',[0.035 0.75 0.7 0.1 ],...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Tooltip', 'Export path',...
                'Callback',@this.edit_dir_export_Callback,...
                'Tag','edit_dir_export');
            %Text Variable to export
            h10 = uicontrol(...
                'Parent',h1,...
                'Units','normalized',...
                'FontUnits',get(0,'defaultuicontrolFontUnits'),...
                'HorizontalAlignment','left',...
                'String','Select variables to export:',...
                'Tooltip', 'Select the variables you want to export',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Style','text',...
                'Position',[0.035 0.63 0.7 0.04 ],...
                'Tag','label_variables'); 
            
            % Table variables to export
            h11 = uitable(...
                'Parent',h1,...
                'Units','normalized',...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontName',matRad_cfg.gui.fontName,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'ColumnWidth',{  20 278 },...
                'Position',[0.035 0.3 0.915 0.3],...
                'ColumnEditable',[true false],...
                'ColumnFormat',{  'logical' 'char' },...
                'RowName',{},...
                'ColumnName',{},...
                'Tag','uitable_variables');
            
            set(h11,'Units','pixels');
            pos = get(h11,'Position');
            set(h11,'ColumnWidth',{20 pos(3)-39});
            
            this.createHandles();
        end
        
    end
    
    
    %CALLBACK METHODS
    methods
        % CALLBACK FOR H2 BUTTON EXPORT
        function this = btn_export_Callback(this,  hObject, event)
            handles = this.handles;
            
            exportDir = get(handles.edit_dir_export,'String');
            
            %Sanity check-
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
               
            try
                dcmExport = matRad_DicomExporter;
                dcmExport.dicomDir = exportDir;
                
                varData = get(handles.uitable_variables,'Data');
                var_selected=false;
                for i= 1:size(varData,1)
                    if varData{i,1} == true
                        var_selected=true;
                        switch varData{i,2}
                            case 'ct'
                                dcmExport.matRad_exportDicomCt();
                            case 'cst'
                                dcmExport.matRad_exportDicomRTStruct();
                            case 'resultGUI'
                                dcmExport.matRad_exportDicomRTDoses();
                        end
                    end
                end
            catch ME
                warning(ME.identifier,'couldn''t export! Reason: %s\n',ME.message)
            end
            
            
            if ~var_selected
                errordlg('No variables selected!');
                return;
            end
            
            
        end
        
        %------------CALLBACK FOR H3 BUTTON CANCEL
        function this = btn_cancel_Callback(this, hObject, event, guidata)
          close;
          % close(handles.figure1);
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
        
    end
end



