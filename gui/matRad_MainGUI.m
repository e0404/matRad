classdef matRad_MainGUI < handle
% matRad_MainGUI   Graphical User Interface configuration class   
% This Class is used to create the main GUI interface where all the widgets
% are called and created.
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
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
        guiHandle
        lastStoragePath = []        
    end
    
    properties (Access = private)
        LogoWidget
        WorkflowWidget
        PlanWidget
        OptimizationWidget
        VisualizationWidget
        ViewerOptionsWidget
        StructureVisibilityWidget
        InfoWidget
        ViewingWidget
        DVHStatsWidget
        eventListeners
        GammaWidget
    end
    
   
    
    methods(Access = protected)
        function this = createMenuBar(this)
           h1 = this.guiHandle;
           folder = fileparts(mfilename('fullpath'));
           loadIcons = load(fullfile(folder,'matRad_iconsGUI.mat'),'icons');

           icons = loadIcons.icons;
            
            h60 = uitoolbar(...
                'Parent',h1,...
                'Tag','uitoolbar1');
            
            h61 = uipushtool(...
                'Parent',h60,...
                'Children',[],...
                'BusyAction','cancel',...
                'Interruptible','off',...
                'Tag','toolbarLoad',...
                'CData',icons{1},...
                'ClickedCallback',@(hObject,eventdata) toolbarLoad_ClickedCallback(this,hObject,eventdata),...
                'Separator','on',...
                'TooltipString','Open File' );
            
            h62 = uipushtool(...
                'Parent',h60,...
                'BusyAction','cancel',...
                'Interruptible','off',...
                'Tag','toolbarSave',...
                'CData',icons{2},...
                'ClickedCallback',@(hObject,eventdata) toolbarSave_ClickedCallback(this,hObject,eventdata),...
                'Separator','on',...
                'TooltipString','Save Figure');
            
            
            h63 = uipushtool(...
                'Parent',h60,...
                'Tag','uipushtool_screenshot',...
                'CData',icons{3},...
                'ClickedCallback',@(hObject,eventdata)uipushtool_screenshot_ClickedCallback(this, hObject, eventdata),...
                'TooltipString','Take a screenshot of the current dose or profile plot' );
            
            h64 = uitoggletool(...
                'Parent',h60,...
                'Tag','toolbarZoomIn',...
                'CData',icons{4},...
                'ClickedCallback',@(hObject, eventdata)toolbarZoomIn_ClickedCallback(this, hObject, eventdata),...
                'Separator','on',...
                'TooltipString','Zoom In');
            
            h65 = uitoggletool(...
                'Parent',h60,...
                'Children',[],...
                'Tag','toolbarZoomOut',...
                'CData',icons{5},...
                'ClickedCallback',@(hObject, eventdata)toolbarZoomOut_ClickedCallback(this, hObject, eventdata),...
                'Separator','on',...
                'TooltipString','Zoom Out');
            
            h66 = uitoggletool(...
                'Parent',h60,...
                'Tag','toolbarPan',...
                'CData',icons{6},...
                'ClickedCallback',@(hObject, eventdata)toolbarPan_ClickedCallback(this, hObject, eventdata),...
                'Separator','on',...
                'TooltipString','Pan' );
            
            h67 = uitoggletool(...
                'Parent',h60,...
                'Tag','toolbarCursor',...
                'CData',icons{7},...
                'ClickedCallback',@(hObject, eventdata)toolbarCursor_ClickedCallback(this, hObject, eventdata),...
                'Separator','on',...
                'TooltipString','Data Cursor' );
            
            h68 = uitoggletool(...
                'Parent',h60,...
                'Tag','toolbarLegend',...
                'CData',icons{8},...
                'ClickedCallback',@(hObject, eventdata)toolbarLegend_ClickedCallback(this, hObject, eventdata),...
                'Separator','on',...
                'TooltipString','Insert Legend');
            
            h69 = uitoggletool(...
                'Parent',h60,...
                'Tag','uitoggletool8',...
                'CData',icons{9},...
                'ClickedCallback',@(hObject, eventdata)uitoggletool8_ClickedCallback(this, hObject, eventdata),...
                'Separator','on',...
                'TooltipString','Insert Colorbar' );
% TODO: future dose comparison tool
%             h70 = uipushtool(...
%                 'Parent',h60,...
%                 'Tag','uipushtool_GammaIndex',...
%                 'CData',icons{10},...
%                 'ClickedCallback',@(hObject,eventdata)gammaIndex_ClickedCallback(this, hObject, eventdata),...
%                 'TooltipString','Calculate Gamma Index' );
            
        end
    end
    
    methods
        function obj = matRad_MainGUI(varargin)
            %Panel for Main Widget 
            
            matRad_cfg = MatRad_Config.instance();
            
            obj.guiHandle = figure(...
                'Units','normalized',...
                'Position',[0.005 0.04 0.99 0.9],... %approximate fullscreen position
                'Visible','on',...
                'Color',matRad_cfg.gui.backgroundColor,...  
                'IntegerHandle','off',...
                'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                'MenuBar','none',...
                'Name','matRad (NOT FOR CLINICAL USE!)',...
                'HandleVisibility','callback',...
                'Tag','figure1',...
                'CloseRequestFcn',@(src,hEvent) figure1_CloseRequestFcn(obj,src,hEvent));
            
            %WindowState not available in all versions
            if isprop(obj.guiHandle,'WindowState')
                set(obj.guiHandle,'WindowState','maximized');
            end
                
            pWorkflow = uipanel(...
                'Parent',obj.guiHandle,...                
                'Title','Workflow',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel4',...
                'Clipping','off',...
                'Position',[0.005 0.8 0.425 0.2],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pPlan = uipanel(...
                'Parent',obj.guiHandle,...
                'Title','Plan',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel1',...
                'Clipping','off',...
                'Position',[0.005 0.55 0.425 0.25],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pOptimization = uipanel(...
                'Parent',obj.guiHandle,...
                'Title',strtrim(strjoin({  'Objectives & constraints'; '                        '; '                        ' })),...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel3',...
                'Clipping','off',...
                'Position',[0.005 0.21 0.425 0.34],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pVisualization = uipanel(...
                'Parent',obj.guiHandle,...
                'Title',strtrim(strjoin({  'Visualization'; '             ' })),...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel2',...
                'Clipping','off',...
                'Position',[0.005 0.005 0.425 0.205],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pLogo = uipanel(...
                'Parent',obj.guiHandle,...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel13',...
                'Clipping','off',...
                'Position',[0.44 0.89 0.55 0.1],...
                'HighLightColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'BorderType','none',...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pViewing = uipanel(...
                'Parent',obj.guiHandle,...
                'Title','Viewing',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel11',...
                'Clipping','off',...
                'Position',[0.435 0.005 0.455 0.885],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pViewerOptions = uipanel(...
                'Parent',obj.guiHandle,...
                'Title','Viewer Options',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel_colormapOptions',...
                'Clipping','off',...
                'Position',[0.895 0.445 0.1 0.445],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pStructureVisibility = uipanel(...
                'Parent',obj.guiHandle,...
                'Title','Structure Visibilty',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel10',...
                'Clipping','off',...
                'Position',[0.895 0.11 0.1 0.33],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
            
            pInfo = uipanel(...
                'Parent',obj.guiHandle,...
                'Title','Info',...
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel12',...
                'Clipping','off',...
                'Position',[0.895 0.005 0.1 0.1],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight);
           
            % Initialize Widgets            
            obj.WorkflowWidget = matRad_WorkflowWidget(pWorkflow);
            obj.PlanWidget = matRad_PlanWidget(pPlan);            
            obj.OptimizationWidget = matRad_OptimizationWidget(pOptimization);
            obj.ViewingWidget = matRad_ViewingWidget(pViewing);
            obj.ViewingWidget.scrollHandle =  obj.guiHandle;            
            obj.VisualizationWidget = matRad_VisualizationWidget(obj.ViewingWidget,pVisualization);        
            obj.ViewerOptionsWidget = matRad_ViewerOptionsWidget(obj.ViewingWidget,pViewerOptions);
            obj.StructureVisibilityWidget = matRad_StructureVisibilityWidget(pStructureVisibility);           
            obj.InfoWidget = matRad_InfoWidget(pInfo); % does not need a listener
            obj.LogoWidget = matRad_LogoWidget(pLogo); % does not need a listener
            
            %Initialize event Listeners
            switch matRad_cfg.env
                case 'MATLAB' 
                    %Event listeners triggered on events
                    obj.eventListeners.workflow = addlistener(obj.WorkflowWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.plan = addlistener(obj.PlanWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.optimization = addlistener(obj.OptimizationWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.viewing = addlistener(obj.ViewingWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.plot = addlistener(obj.ViewingWidget,'plotUpdated',@(src,hEvent) updateButtons(obj));
                    obj.eventListeners.visualization = addlistener(obj.VisualizationWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.viewerOptions = addlistener(obj.ViewerOptionsWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.structureVisibility = addlistener(obj.StructureVisibilityWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    
                    % only available in MATLAB
                    obj.ViewingWidget.dcmHandle = datacursormode(obj.guiHandle);
                    obj.ViewingWidget.panHandle = pan(obj.guiHandle);
                    obj.ViewingWidget.zoomHandle = zoom(obj.guiHandle);
                    
                case 'OCTAVE'
                    % addlistener is not yet available in octave 
                    obj.eventListeners.workflow = matRad_addListenerOctave(obj.WorkflowWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.plan = matRad_addListenerOctave(obj.PlanWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.optimization = matRad_addListenerOctave(obj.OptimizationWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.viewing = matRad_addListenerOctave(obj.ViewingWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.plot = matRad_addListenerOctave(obj.ViewingWidget,'plotUpdated',@(src,hEvent) updateButtons(obj));
                    obj.eventListeners.visualization = matRad_addListenerOctave(obj.VisualizationWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.viewerOptions = matRad_addListenerOctave(obj.ViewerOptionsWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    obj.eventListeners.structureVisibility = matRad_addListenerOctave(obj.StructureVisibilityWidget,'workspaceChanged',@(src,hEvent) updateWidgets(obj,hEvent));
                    
            end
            
            obj.createMenuBar();
            
            % update button states
            obj.updateButtons();
            
            try
                % change color of toobar the first time GUI is started
                hToolbar = findall(obj.guiHandle,'tag','uitoolbar1');
                jToolbar = get(get(hToolbar,'JavaContainer'),'ComponentPeer');
                jToolbar.setBorderPainted(false);
                color = java.awt.Color.gray;
                % Remove the toolbar border, to blend into figure contents
                jToolbar.setBackground(color);
                % Remove the separator line between toolbar and contents
                warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
                jFrame = get(handle(obj.guiHandle),'JavaFrame');
                jFrame.showTopSeparator(false);
                jtbc = jToolbar.getComponents;
                for idx=1:length(jtbc)
                    jtbc(idx).setOpaque(false);
                    jtbc(idx).setBackground(color);
                    for childIdx = 1 : length(jtbc(idx).getComponents)
                        jtbc(idx).getComponent(childIdx-1).setBackground(color);
                    end
                end
            catch
                matRad_cfg.dispDeprecationWarning('Java properties couldn''t be set');
            end
        end
        
        function this = updateWidgets(this,evt)
           matRad_cfg = MatRad_Config.instance();
           if matRad_cfg.logLevel > 3 
               if nargin < 2 || ~isa(evt,'matRad_WorkspaceChangedEvent')
                   changed = {};
               else
                   changed = evt.changedVariables;
               end

               matRad_cfg.dispDebug('GUI Workspace Changed at %s. ',datestr(now,'HH:MM:SS.FFF')); 
               if ~isempty(changed)
                   matRad_cfg.dispDebug('Specific Variables: %s.\n',strjoin(changed,'|'));
               else
                   matRad_cfg.dispDebug('\n');
               end
           end
           
          
           %%Debug
           this.PlanWidget.update(evt);
           this.WorkflowWidget.update(evt);
           this.OptimizationWidget.update(evt);
           this.ViewingWidget.lockUpdate = 0;
           this.ViewingWidget.update(evt);
           this.StructureVisibilityWidget.update(evt);
           %this.ViewerOptionsWidget.update();
           %this.VisualizationWidget.update();
           
        end
        
        function this = updateButtons(this)
           %disp(['Plot Changed ' datestr(now,'HH:MM:SS.FFF')]); %Debug
           % update the visualization and viewer options widgets
           
           matRad_cfg = MatRad_Config.instance();

           if strcmp(matRad_cfg.env,'OCTAVE')
              return
           end
         
           set(findobj(this.guiHandle,'tag','toolbarPan'),'State',get(this.ViewingWidget.panHandle,'Enable'));
           set(findobj(this.guiHandle,'tag','toolbarCursor'),'State',get(this.ViewingWidget.dcmHandle,'Enable'));
             
           if this.ViewingWidget.plotColorBar  && ~isempty(this.ViewingWidget.cBarHandle) && isvalid(this.ViewingWidget.cBarHandle)
               set(findobj(this.guiHandle,'tag','uitoggletool8'),'State','on')
           else
               set(findobj(this.guiHandle,'tag','uitoggletool8'),'State','off')
           end
           if this.ViewingWidget.plotLegend && ~isempty(this.ViewingWidget.legendHandle) && isvalid(this.ViewingWidget.legendHandle)
               set(findobj(this.guiHandle,'tag','toolbarLegend'),'State',get(this.ViewingWidget.legendHandle,'visible'));
           else
               set(findobj(this.guiHandle,'tag','toolbarLegend'),'State','off');
           end
           if strcmp(get(this.ViewingWidget.zoomHandle,'Enable'),'on')
               if strcmp(get(this.ViewingWidget.zoomHandle,'Direction'),'in')
                   set(findobj(this.guiHandle,'tag','toolbarZoomOut'),'State','off');
                   set(findobj(this.guiHandle,'tag','toolbarZoomIn'),'State','on');
               else
                   set(findobj(this.guiHandle,'tag','toolbarZoomOut'),'State','on');
                   set(findobj(this.guiHandle,'tag','toolbarZoomIn'),'State','off');
               end
           else
               set(findobj(this.guiHandle,'tag','toolbarZoomOut'),'State','off');
               set(findobj(this.guiHandle,'tag','toolbarZoomIn'),'State','off');
           end
           
           this.ViewerOptionsWidget.update();
           this.VisualizationWidget.update();
           
        end
        
        %% Callbacks
        % toolbar load button
        function toolbarLoad_ClickedCallback(this,hObject, eventdata)
            this.WorkflowWidget.btnLoadMat_Callback(hObject, eventdata);
        end
        % toolbar save button
        function toolbarSave_ClickedCallback(this,hObject, eventdata)
            %handles=this.handles;

            
            answer = questdlg('Do you wish to save the full workspace or only matRad variables?','Save','Full workspace', 'matRad variables', 'matRad variables');
            
            
            switch answer
                case 'Full workspace'
                    uisave;
                case 'matRad variables'
                    variables = {'cst','ct','pln','stf','dij','resultGUI'};
                    vExists=false(size(variables));
                    for i= 1:numel(variables)
                        var= char(variables(i));
                        vExists(i) = evalin('base',['exist(''' var ''',''var'')']);
                        if vExists(i)
                            eval([var '=evalin(''base'',''' var ''');'])
                        end 
                    end
                    uisave(variables(vExists));
            end

        end
        %Toolbar screenshot of the Viewing widget
        function uipushtool_screenshot_ClickedCallback(this,hObject, eventdata)
            % hObject    handle to uipushtool_screenshot (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            
            tmpFig = figure('position',[100 100 700 600],'Visible','off','name','Current View');
            cBarHandle = this.ViewingWidget.cBarHandle; %findobj(handles.figure1,'Type','colorbar');
            if ~isempty(cBarHandle)
                new_handle = copyobj([this.ViewingWidget.handles.axesFig cBarHandle],tmpFig);                
            else
                new_handle = copyobj(this.ViewingWidget.handles.axesFig,tmpFig);                
            end
            
            oldPos = get(this.ViewingWidget.handles.axesFig,'Position');
            set(new_handle(1),'units','normalized', 'Position',oldPos);
            
            if exist(this.lastStoragePath,'dir') ~= 7 
                this.lastStoragePath = [];
            end
            
            [filename, pathname] = uiputfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'; '*.fig','MATLAB figure file'},'Save current view',[this.lastStoragePath 'screenshot.png']);
            
            this.lastStoragePath = pathname;
            
            if ~isequal(filename,0) && ~isequal(pathname,0)
                set(gcf, 'pointer', 'watch');
                saveas(tmpFig,fullfile(pathname,filename));
                set(gcf, 'pointer', 'arrow');
                close(tmpFig);
                uiwait(msgbox('Current view has been succesfully saved!'));
            else
                uiwait(msgbox('Aborted saving, showing figure instead!'));
                set(tmpFig,'Visible','on');
            end
            
        end
       
        % Toggle color bar callback
        function uitoggletool8_ClickedCallback(this,hObject, eventdata)
            % hObject    handle to uitoggletool8 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            this.ViewingWidget.plotColorBar = strcmp(get(hObject,'State'),'on');

        end
        %Toggle Legend Callback
        function toolbarLegend_ClickedCallback(this,hObject, eventdata)
           
            this.ViewingWidget.plotLegend = strcmp(get(hObject,'State'),'on');
        end
        %Toggle Zoom in Callback
        function toolbarZoomIn_ClickedCallback(this,hObject, eventdata)
            set(this.ViewingWidget.zoomHandle,'Enable',char(get(hObject,'State')));
            set(this.ViewingWidget.zoomHandle,'Direction','in');
            set(findobj('tag','toolbarZoomOut'),'State','off');
        end
        %Toggle Zoom out Callback
        function toolbarZoomOut_ClickedCallback(this,hObject, eventdata)
            set(this.ViewingWidget.zoomHandle,'Enable',char(get(hObject,'State')));
            set(this.ViewingWidget.zoomHandle,'Direction','out');
            set(findobj('tag','toolbarZoomIn'),'State','off');
        end
        %Toggle Pan Callback
        function toolbarPan_ClickedCallback(this,hObject, eventdata)
           set(this.ViewingWidget.panHandle,'Enable',char(get(hObject,'State')));
        end
        %Toggle cursor Callback
        function toolbarCursor_ClickedCallback(this,hObject, eventdata)
           set(this.ViewingWidget.dcmHandle,'Enable',get(hObject,'State'));
        end        

        %Toggle cursor Callback
        function gammaIndex_ClickedCallback(this,hObject, eventdata)
           % this.GammaWidget = matRad_GammaWidget();
        end     

        
        % button: close
        function figure1_CloseRequestFcn(this,hObject, ~)
            matRad_cfg = MatRad_Config.instance();
            set(0,'DefaultUicontrolBackgroundColor',matRad_cfg.gui.backgroundColor);
            selection = questdlg('Do you really want to close matRad?',...
                'Close matRad',...
                'Yes','No','Yes');

            switch selection
                case 'Yes'
                    delete(hObject);
                    delete(this);
                case 'No'
                    return
            end
        end
        
    end
    
    
end

