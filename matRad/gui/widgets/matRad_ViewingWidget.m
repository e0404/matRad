classdef matRad_ViewingWidget < matRad_Widget
    % matRad_ViewingWidget class to generate GUI widget to display plan
    % dose distributions and ct
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
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        plane = 3;
        slice = 1;
        selectedBeam = 1;
        numOfBeams=1;
        profileOffset=0;
        OffsetSliderStep;
        OffsetMinMax;
        typeOfPlot = 1;
        colorData;
        doseColorMap = 'jet';
        ctColorMap = 'bone';
        cMapSize = 64;
        ctScen = 1;
        plotCT = true;
        plotContour = true;
        plotIsoCenter = true;
        plotPlan = false;
        plotDose = true;
        plotIsoDoseLines = true;
        plotIsoDoseLinesLabels = false;
        plotLegend = false;
        plotColorBar = true;
        ProfileType = 'lateral';
        SelectedDisplayOption = 'physicalDose';
        SelectedDisplayAllOptions = '';
        CutOffLevel = 0.01;
        dispWindow = cell(3,2);
        doseOpacity = 0.6;
        IsoDose_Levels= [];
        NewIsoDoseFlag = true;
        cBarHandle;
        dcmHandle;
        panHandle;
        zoomHandle;
        legendHandle;
        scrollHandle;
        lockColorSettings = false;
        %plotlegend=false;
        evt;
    end
    
    properties (SetAccess=private)
        
        IsoDose_Contours; %only updated from within this class
        VOIPlotFlag;
        DispInfo;
        AxesHandlesVOI;        
        cst;
        vIsoCenter;
        sliceContourLegend;
    end

    properties (SetAccess = protected)
        initialized = false;
    end
    
    events
        plotUpdated
    end
    
    methods
        function this = matRad_ViewingWidget(handleParent)
            matRad_cfg = MatRad_Config.instance();

            if nargin < 1                
                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.3 0.2 0.4 0.6],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,...  'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Viewing',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1');
                
            end
            
            this = this@matRad_Widget(handleParent);
            
            matRad_cfg = MatRad_Config.instance();
            
            if nargin < 1
                % create the handle objects if there's no parent
                this.scrollHandle = this.widgetHandle;
                
                % only available in MATLAB
                if matRad_cfg.isMatlab
                    this.dcmHandle = datacursormode(this.widgetHandle);
                    this.panHandle = pan(this.widgetHandle);
                    this.zoomHandle = zoom(this.widgetHandle);
                end
            end
            this.initialize();
        end

        function initialize(this)
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
                        
        function notifyPlotUpdated(obj)
            % handle environment
            matRad_cfg = MatRad_Config.instance();
            switch matRad_cfg.env
                case 'MATLAB'
                    notify(obj, 'plotUpdated');
                case 'OCTAVE'
                    matRad_notifyOctave(obj, 'plotUpdated');
            end
            
        end
        
        %% SET FUNCTIONS
        function set.plane(this,value)
            this.plane=value;
            this.update();
        end
        
        function set.slice(this,value)
            % project to allowed set (between min and max value)
            newSlice = max(value,1);
            if evalin('base','exist(''ct'')') 
                ct=evalin('base','ct');
                newSlice = min(newSlice,ct.cubeDim(this.plane));
            else
                newSlice=1;
            end
            
            this.slice=newSlice;
            this.update();
        end
                       
        function set.selectedBeam(this,value)
            this.selectedBeam=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.numOfBeams(this,value)
            this.numOfBeams=value;
            this.update();
        end
        
        function set.profileOffset(this,value)
            this.profileOffset=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.OffsetSliderStep(this,value)
            this.OffsetSliderStep=value;
            this.update();
        end
        
        function set.OffsetMinMax(this,value)
            this.OffsetMinMax=value;
            this.update();
        end
        
        
        function set.typeOfPlot(this,value)
            this.typeOfPlot=value;            
            evt = matRad_WorkspaceChangedEvent('image_display');
            cla(this.handles.axesFig,'reset');
            this.update(evt);
        end
        
        function set.colorData(this,value)
            this.colorData=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.doseColorMap(this,value)
            this.doseColorMap=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.ctColorMap(this,value)
            this.ctColorMap=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.cMapSize(this,value)
            this.cMapSize=value;
            this.update();
        end

        function set.ctScen(this,value)
            this.ctScen=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotCT(this,value)
            this.plotCT=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotContour(this,value)
            this.plotContour=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotIsoCenter(this,value)
            this.plotIsoCenter=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotPlan(this,value)
            this.plotPlan=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotDose(this,value)
            this.plotDose=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotIsoDoseLines(this,value)
            this.plotIsoDoseLines=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotIsoDoseLinesLabels(this,value)
            this.plotIsoDoseLinesLabels=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotLegend(this,value)
            this.plotLegend=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.plotColorBar(this,value)
            this.plotColorBar=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
                
        function set.ProfileType(this,value)
            this.ProfileType=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.SelectedDisplayOption(this,value)
            this.SelectedDisplayOption=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.SelectedDisplayAllOptions(this,value)
            this.SelectedDisplayAllOptions=value;
            this.update();
        end
        
        
        function set.CutOffLevel(this,value)
            this.CutOffLevel=value;
            this.update();
        end
        
        function set.dispWindow(this,value)
            this.dispWindow=value;
            evt = matRad_WorkspaceChangedEvent('viewer_options');
            this.update(evt);
        end
        
        function set.doseOpacity(this,value)
            this.doseOpacity=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.IsoDose_Levels(this,value)
            this.IsoDose_Levels=value;
            evt = matRad_WorkspaceChangedEvent('viewer_options');
            this.update(evt);
        end
        
        function set.IsoDose_Contours(this,value)
            this.IsoDose_Contours=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.NewIsoDoseFlag(this,value)
            this.NewIsoDoseFlag=value;
            evt = matRad_WorkspaceChangedEvent('image_display');
            this.update(evt);
        end
        
        function set.dcmHandle(this,value)
            this.dcmHandle=value;
            set(this.dcmHandle,'DisplayStyle','window');
            %Add the callback for the datacursor display
            set(this.dcmHandle,'UpdateFcn',@(hObject, eventdata)dataCursorUpdateFunction(this, hObject, eventdata));
        end
        
        function set.scrollHandle(this,value)
            this.scrollHandle=value;
            % set callback for scroll wheel function
            set(this.scrollHandle,'WindowScrollWheelFcn',@(src,event)matRadScrollWheelFcn(this,src,event));
        end
    end
 %%
    methods(Access = protected)
        function this = createLayout(this)
            %Viewer Widget 
            h88 = this.widgetHandle;

            matRad_cfg = MatRad_Config.instance();
                        
            h89 = axes(...
                'Parent',h88,...
                'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
                'XTickLabel',{  '0'; '0.1'; '0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1' },...
                'XColor',matRad_cfg.gui.textColor,...
                'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
                'YTickLabel',{  '0'; '0.1'; '0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1' },...
                'YColor',matRad_cfg.gui.textColor,...
                'ZColor',matRad_cfg.gui.textColor,...
                'Units','normalized',...
                'Position',[0.0718390804597701 0.0654391371340524 0.902298850574712 0.899121725731895],...
                'Color',matRad_cfg.gui.elementColor,...
                'Box','on',...
                'BoxStyle','full',...
                'Tag','axesFig'); 
                 
            %Title
            h90 = get(h89,'title');
            
            set(h90, ...
                'Parent',h89,...
                'Visible','on',...
                'Units','data',...
                'Position',[0.500000554441759 1.00453467465753 0.5],...
                'PositionMode','auto', ...
                'Margin',2,...
                'Clipping','off',...
                'FontUnits','points',...
                'Color',[0 0 0],...
                'BackgroundColor','none',...
                'EdgeColor','none',...
                'Interpreter','tex',...
                'Rotation',0,...
                'RotationMode','auto',...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontAngle','normal',...
                'FontWeight','normal',...
                'HorizontalAlignment','center',...
                'HorizontalAlignmentMode','auto',...
                'VerticalAlignment','bottom',...
                'VerticalAlignmentMode','auto', ...
                'LineStyle','-',...
                'LineWidth',0.5,...
                'HandleVisibility','off',...
                'BusyAction','queue',...
                'Interruptible','on',...
                'HitTest','on');

            %X Label
            h91 = get(h89,'xlabel');
            
            set(h91,...
                'Parent',h89,...
                'Units','data',...
                'FontUnits','points',...
                'Color',matRad_cfg.gui.textColor,...
                'Position',[0.500000476837158 -0.0373767115122652 0],...
                'PositionMode','auto',...
                'Interpreter','tex',...
                'Rotation',0,...
                'RotationMode','auto',...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontAngle','normal',...
                'FontWeight','normal',...
                'HorizontalAlignment','center',...
                'HorizontalAlignmentMode','auto',...
                'VerticalAlignment','top',...
                'VerticalAlignmentMode','auto',...
                'EdgeColor','none',...
                'LineStyle','-',...
                'LineWidth',0.5,...
                'BackgroundColor','none',...
                'Margin',3,...
                'Clipping','off',...
                'Visible','on',...
                'HandleVisibility','off',...
                'BusyAction','queue',...
                'Interruptible','on',...
                'HitTest','on');
            
            %Y label
            h92 = get(h89,'ylabel');
            
            set(h92,...
                'Parent',h89,...
                'Units','data',...
                'FontUnits','points',...
                'Color',matRad_cfg.gui.textColor,...
                'Position',[-0.0474647368237942 0.500000476837158 0],...
                'PositionMode','auto',...
                'Interpreter','tex',...
                'Rotation',90,...
                'RotationMode','auto',...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontAngle','normal',...
                'FontWeight','normal',...
                'HorizontalAlignment','center',...
                'HorizontalAlignmentMode','auto',...
                'VerticalAlignment','bottom',...
                'VerticalAlignmentMode','auto',...
                'EdgeColor','none',...
                'LineStyle','-',...
                'LineWidth',0.5,...
                'BackgroundColor','none',...
                'Margin',3,...
                'Clipping','off',...
                'Visible','on',...
                'HandleVisibility','off',...
                'BusyAction','queue',...
                'Interruptible','on',...
                'HitTest','on');
            
            %Z label
            h93 = get(h89,'zlabel');
            
            set(h93,...
                'Parent',h89,...
                'Units','data',...
                'FontUnits','points',...
                'Color',matRad_cfg.gui.textColor,...
                'Position',[0 0 0],...
                'PositionMode','auto',...
                'Interpreter','tex',...
                'Rotation',0,...
                'RotationMode','auto',...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontAngle','normal',...
                'FontWeight','normal',...
                'HorizontalAlignment','left',...
                'HorizontalAlignmentMode','auto',...
                'VerticalAlignment','middle',...
                'VerticalAlignmentMode','auto',...
                'EdgeColor','none',...
                'LineStyle','-',...
                'LineWidth',0.5,...
                'BackgroundColor','none',...
                'Margin',3,...
                'Clipping','off',...
                'Visible','off',...
                'HandleVisibility','off',...
                'BusyAction','queue',...
                'Interruptible','on',...
                'HitTest','on');           
            
            this.createHandles();
            
        end
    
        function this=doUpdate(this,evt)
            if ~this.initialized
                this.initValues();
                this.initialized = true;
            end

            if ~this.updateLock
            
                %doUpdate = false;
                if nargin == 2
                    this.evt = evt;
                    %At pln changes and at cst/cst (for Isocenter and new settings) 
                    %we need to update
                    if this.checkUpdateNecessary({'ct','cst'},evt)
                        this.initValues();
                    end
                    this.updateValues();
                    %doUpdate = this.checkUpdateNecessary({'pln_display','ct','cst','resultGUI','image_display'},evt);
                    if  this.checkUpdateNecessary({'resultGUI','image_display','viewer_options','pln_angles'},evt)
                        this.updateIsoDoseLineCache();
                        this.UpdatePlot();
                    end                 
                else
                    this.updateValues();
                    this.UpdatePlot();
                end
                            
                this.evt =[];
            end
            
        end
    end
    
    methods
        function UpdatePlot(this)
            if this.updateLock
                return
            end

            matRad_cfg = MatRad_Config.instance();
            
            handles = this.handles;
            
            %profile on;
            axes(handles.axesFig);
            
            % this is necessary to prevent multiple callbacks of update plot drawing on
            % top of each other in matlab <2014
            drawnow;
            
            defaultFontSize     = matRad_cfg.gui.fontSize;
            currAxes            = axis(handles.axesFig);
            axesHandlesVOI      = cell(0);
            
            axesHandlesCT_Dose  = cell(0);
            axesHandlesIsoDose  = cell(0);
            
            if evalin('base','exist(''ct'')')
                 ct = evalin('base','ct');
            else
                 cla(handles.axesFig);
                 return
            end

            if evalin('base','exist(''pln'')')
                pln = evalin('base','pln');
            else
                pln = [];
            end
            
            %% If resultGUI exists, then an optimization has been performed
            if evalin('base','exist(''resultGUI'')') 
                result = evalin('base','resultGUI');
            end

            %% set and get required variables
            %plane = get(handles.popupPlane,'Value');
            %slice = round(get(handles.sliderSlice,'Value'));
            hold(handles.axesFig,'on');
            if this.typeOfPlot==1 %get(handles.popupTypeOfPlot,'Value')==1
                set(handles.axesFig,'YDir','Reverse');
            end
                 
            selectIx = this.colorData; %get(handles.popupmenu_chooseColorData,'Value');
            
            cla(handles.axesFig);
            %% plot ct - if a ct cube is available and type of plot is set to 1 and not 2; 1 indicate cube plotting and 2 profile plotting
            if ~isempty(ct) && this.typeOfPlot==1
                
                if selectIx == 2
                    ctIx = 1;
                else
                    ctIx = selectIx;
                end
                if isfield(ct, 'cubeHU')
                    plotCtCube = ct.cubeHU;
                else
                    plotCtCube = ct.cube;
                end
                ctMap = matRad_getColormap(this.ctColorMap,this.cMapSize);
                
%                 if isempty(this.dispWindow{ctIx,2})
%                     this.dispWindow{ctIx,2} = [min(reshape([ct.cubeHU{:}],[],1)) max(reshape([ct.cubeHU{:}],[],1))];
%                 end
                
                if this.plotCT %get(handles.radiobtnCT,'Value')
                    [axesHandlesCT_Dose{end+1},~,~] = matRad_plotCtSlice(handles.axesFig,plotCtCube,this.ctScen,this.plane,this.slice,ctMap,this.dispWindow{ctIx,1});
                    
                    % plot colorbar? If 1 the user asked for the CT.
                    % not available in octave 
                    if strcmp(matRad_cfg.env,'MATLAB') && ~isempty(this.colorData) && this.colorData == 1
                        %Plot the colorbar
                        this.cBarHandle = matRad_plotColorbar(handles.axesFig,ctMap,this.dispWindow{ctIx,1},'FontSize',defaultFontSize,'Color',matRad_cfg.gui.textColor);
                        
                        if this.plotColorBar
                            set(this.cBarHandle,'Visible','on')
                        else
                            set(this.cBarHandle,'Visible','off')
                        end
                        
                        %adjust lables
                        if isfield(ct,'cubeHU')
                            set(get(this.cBarHandle,'ylabel'),'String', 'Hounsfield Units','fontsize',defaultFontSize);
                        else
                            set(get(this.cBarHandle,'ylabel'),'String', 'Electron Density','fontsize',defaultFontSize);
                        end
                        % do not interprete as tex syntax
                        set(get(this.cBarHandle,'ylabel'),'interpreter','none');
                    end
                end
            end
            
            %% plot dose cube
            if this.typeOfPlot== 1  && exist('result','var') % handles.State >= 1 && 
                doseMap = matRad_getColormap(this.doseColorMap,this.cMapSize);
                doseIx  = 2;
                
                dose = result.(this.SelectedDisplayOption);
                
                % dose colorwash
                if ~isempty(dose) && ~isvector(dose)
                  
                    if this.plotDose 
                        [doseHandle,~,~] = matRad_plotDoseSlice(handles.axesFig,dose,this.plane,this.slice,this.CutOffLevel,this.doseOpacity,doseMap,this.dispWindow{doseIx,1});
                        axesHandlesCT_Dose{end+1}         = doseHandle;
                    end
                    
                    % plot colorbar
                    if matRad_cfg.isMatlab && ~isempty(this.colorData) && this.colorData > 1 
                        %Plot the colorbar
                        this.cBarHandle = matRad_plotColorbar(handles.axesFig,doseMap,this.dispWindow{selectIx,1},'fontsize',defaultFontSize,'Color',matRad_cfg.gui.textColor);
                        
                        if this.plotColorBar
                            set(this.cBarHandle,'Visible','on')
                        else
                            set(this.cBarHandle,'Visible','off')
                        end
                        %adjust lables
                        Idx = find(strcmp(this.SelectedDisplayOption,this.DispInfo(:,1)));
                        set(get(this.cBarHandle,'ylabel'),'String', [this.DispInfo{Idx,1} ' ' this.DispInfo{Idx,3} ],'fontsize',defaultFontSize);
                        % do not interprete as tex syntax
                        set(get(this.cBarHandle,'ylabel'),'interpreter','none');
                    end
                end
                
                
                %% plot iso dose lines
                if this.plotIsoDoseLines 
                    plotLabels = this.plotIsoDoseLinesLabels; 
                    axesHandlesIsoDose = matRad_plotIsoDoseLines(handles.axesFig,dose,this.IsoDose_Contours,this.IsoDose_Levels,plotLabels,this.plane,this.slice,doseMap,this.dispWindow{doseIx,1},'LineWidth',1.5);
                end
            end
           
            %% plot VOIs
            if this.plotContour && this.typeOfPlot==1 && exist('ct','var') %&& get(handles.radiobtnContour,'Value') && handles.State>0
                [AxVOI, this.sliceContourLegend] = matRad_plotVoiContourSlice(handles.axesFig,this.cst,ct,this.ctScen,this.VOIPlotFlag,this.plane,this.slice,[],'LineWidth',2);
                axesHandlesVOI = [axesHandlesVOI AxVOI];
            end
            this.AxesHandlesVOI=axesHandlesVOI;
            
            %% Set axis labels and plot iso center
            matRad_plotAxisLabels(handles.axesFig,ct,this.plane,this.slice,defaultFontSize);
            set(get(handles.axesFig,'Title'),'Color',matRad_cfg.gui.textColor);
                        
            if this.plotIsoCenter && this.typeOfPlot == 1 && ~isempty(pln) %get(handles.radioBtnIsoCenter,'Value') == 1 
                hIsoCenterCross = matRad_plotIsoCenterMarker(handles.axesFig,pln,ct,this.plane,this.slice);
            end
            
            if this.plotPlan && ~isempty(pln) %get(handles.radiobtnPlan,'value') == 1
                matRad_plotProjectedGantryAngles(handles.axesFig,pln,ct,this.plane);
            end
            
            %set axis ratio
            ratios = [1/ct.resolution.x 1/ct.resolution.y 1/ct.resolution.z];
            set(handles.axesFig,'DataAspectRatioMode','manual');
            if this.plane == 1
                res = [ratios(3) ratios(1)]./max([ratios(3) ratios(1)]);
                set(handles.axesFig,'DataAspectRatio',[res 1])
            elseif this.plane == 2 % sagittal plane
                res = [ratios(3) ratios(2)]./max([ratios(3) ratios(2)]);
                set(handles.axesFig,'DataAspectRatio',[res 1])
            elseif  this.plane == 3 % Axial plane
                res = [ratios(1) ratios(2)]./max([ratios(1) ratios(2)]);
                set(handles.axesFig,'DataAspectRatio',[res 1])
            end

            axis(handles.axesFig,'tight');
            
            
            %% profile plot
            if this.typeOfPlot == 2 && exist('result','var')
                % set SAD
                fileName = [pln.radiationMode '_' pln.machine];
                try
                    load(fileName);
                    SAD = machine.meta.SAD;
                catch
                    this.showError(['Could not find the following machine file: ' fileName ]);
                end
                
                % clear view and initialize some values
                cla(handles.axesFig,'reset')
                set(handles.axesFig,'YDir','normal','Color',matRad_cfg.gui.elementColor,'XColor',matRad_cfg.gui.textColor);
                hold(handles.axesFig,'on');
                tmpColor = rgb2hsv(matRad_cfg.gui.elementColor);
                if tmpColor(3) > 0.5
                    tmpColor = 'black';
                else
                    tmpColor = 'white';
                end
                ylabel(['{\color{' tmpColor '}dose [Gy]}']);
                cColor={tmpColor,'green','magenta','cyan','yellow','red','blue'};
                
                % Rotate the system into the beam.
                % passive rotation & row vector multiplication & inverted rotation requires triple matrix transpose
                rotMat_system_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles(this.selectedBeam),pln.propStf.couchAngles(this.selectedBeam)));
                
                if strcmp(this.ProfileType,'longitudinal')
                    sourcePointBEV = [this.profileOffset -SAD   0];
                    targetPointBEV = [this.profileOffset  SAD   0];
                elseif strcmp(this.ProfileType,'lateral')
                    sourcePointBEV = [-SAD this.profileOffset   0];
                    targetPointBEV = [ SAD this.profileOffset   0];
                end
                
                rotSourcePointBEV = sourcePointBEV * rotMat_system_T;
                rotTargetPointBEV = targetPointBEV * rotMat_system_T;
                
                % perform raytracing on the central axis of the selected beam, use unit
                % electron density for plotting against the geometrical depth
                cubeIsoCenter = matRad_world2cubeCoords(pln.propStf.isoCenter(this.selectedBeam,:),ct);
                [~,l,rho,~,ix] = matRad_siddonRayTracer(cubeIsoCenter,ct.resolution,rotSourcePointBEV,rotTargetPointBEV,{0*ct.cubeHU{1}+1});
                d = [0 l .* rho{1}];
                % Calculate accumulated d sum.
                vX = cumsum(d(1:end-1));
               
                % plot physical dose
                %Content =this.SelectedDisplayOption; %get(this.popupDisplayOption,'String');
                SelectedCube = this.SelectedDisplayOption; %Content{get(this.popupDisplayOption,'Value')};
                if sum(strcmp(SelectedCube,{'physicalDose','effect','RBExDose','alpha','beta','RBE','BED'})) > 0
                    Suffix = '';
                else
                    Idx    = find(SelectedCube == '_');
                    Suffix = SelectedCube(Idx:end);
                end
                
                mPhysDose = result.(['physicalDose' Suffix]);
                PlotHandles{1} = plot(handles.axesFig,vX,mPhysDose(ix),'color',cColor{1,1},'LineWidth',3); hold(handles.axesFig,'on');
                PlotHandles{1,2} ='physicalDose';
                ylabel(handles.axesFig,'dose in [Gy]');
                set(handles.axesFig,'FontSize',defaultFontSize);
                 
                % plot counter
                Cnt=2;

                rightAx = [];
                
                if isfield(result,['RBE' Suffix])
                    
                    %disbale specific plots
                    %this.DispInfo{6,2}=0;
                    %this.DispInfo{5,2}=0;
                    %this.DispInfo{2,2}=0;
                    
                    % generate two lines for ylabel
                    StringYLabel1 = ['\fontsize{8}{\color{red}RBE x dose [Gy(RBE)] \color{' tmpColor '}dose [Gy] '];
                    StringYLabel2 = '';
                    for i=1:1:size(this.DispInfo,1)
                        if this.DispInfo{i,2} && sum(strcmp(this.DispInfo{i,1},{['effect' Suffix],['alpha' Suffix],['beta' Suffix]})) > 0
                            %physicalDose is already plotted and RBExDose vs RBE is plotted later with plotyy
                            if ~strcmp(this.DispInfo{i,1},['RBExDose' Suffix]) &&...
                                    ~strcmp(this.DispInfo{i,1},['RBE' Suffix]) && ...
                                    ~strcmp(this.DispInfo{i,1},['physicalDose' Suffix])
                                
                                mCube = result.([this.DispInfo{i,1}]);
                                PlotHandles{Cnt,1} = plot(handles.axesFig,vX,mCube(ix),'color',cColor{1,Cnt},'LineWidth',3); hold(handles.axesFig,'on');
                                PlotHandles{Cnt,2} = this.DispInfo{i,1};
                                StringYLabel2 = [StringYLabel2  ' \color{'  cColor{1,Cnt} '}' this.DispInfo{i,1} ' ['  this.DispInfo{i,3} ']'];
                                Cnt = Cnt+1;
                            end
                        end
                    end
                    StringYLabel2 = [StringYLabel2 '}'];
                    % always plot RBExDose against RBE
                    mRBExDose = result.(['RBExDose' Suffix]);
                    vBED = mRBExDose(ix);
                    mRBE = result.(['RBE' Suffix]);
                    vRBE = mRBE(ix);
                    
                    % plot biological dose against RBE
                    [ax, PlotHandles{Cnt,1}, PlotHandles{Cnt+1,1}]=plotyy(handles.axesFig,vX,vBED,vX,vRBE,'plot'); 
                    hold(ax(2),'on');
                    PlotHandles{Cnt,2}='RBExDose';
                    PlotHandles{Cnt+1,2}='RBE';                    
                    
                    % set plotyy properties
                    set(get(ax(2),'Ylabel'),'String','RBE [1]','FontSize',defaultFontSize);

                    ylabel({StringYLabel1;StringYLabel2})
                    set(PlotHandles{Cnt,1},'Linewidth',4,'color','r');
                    set(PlotHandles{Cnt+1,1},'Linewidth',3,'color','b');
                    set(ax(1),'ycolor','r')
                    set(ax(2),'ycolor','b')
                    set(ax,'FontSize',defaultFontSize);
                    Cnt=Cnt+2;
                    rightAx = ax(2);
                end
                
                % asses target coordinates
                tmpPrior = intmax;
                tmpSize = 0;
                for i=1:size(this.cst,1)
                    if strcmp(this.cst{i,3},'TARGET') && tmpPrior >= this.cst{i,5}.Priority && tmpSize<numel(this.cst{i,4}{1})
                        linIdxTarget = unique(this.cst{i,4}{1});
                        tmpPrior=this.cst{i,5}.Priority;
                        tmpSize=numel(this.cst{i,4}{1});
                        VOI = this.cst{i,2};
                    end
                end
                
                str = sprintf('profile plot - central axis of %d beam, gantry angle %d°, couch angle %d°',...
                    this.selectedBeam ,pln.propStf.gantryAngles(this.selectedBeam),pln.propStf.couchAngles(this.selectedBeam));
                h_title = title(handles.axesFig,str,'FontSize',defaultFontSize,'Color',matRad_cfg.gui.highlightColor);
                pos = get(h_title,'Position');
                set(h_title,'Position',[pos(1)-40 pos(2) pos(3)])
                
                % plot target boundaries
                mTargetCube = zeros(ct.cubeDim);
                mTargetCube(linIdxTarget) = 1;
                vProfile = mTargetCube(ix);
                WEPL_Target_Entry = vX(find(vProfile,1,'first'));
                WEPL_Target_Exit  = vX(find(vProfile,1,'last'));
                PlotHandles{Cnt,2} =[VOI ' boundary'];
                
                if ~isempty(WEPL_Target_Entry) && ~isempty(WEPL_Target_Exit)
                    hold(handles.axesFig,'on');
                    PlotHandles{Cnt,1} = ...
                        plot([WEPL_Target_Entry WEPL_Target_Entry],get(handles.axesFig,'YLim'),'--','Linewidth',3,'color',matRad_cfg.gui.highlightColor);hold(handles.axesFig,'on');
                    plot([WEPL_Target_Exit WEPL_Target_Exit],get(handles.axesFig,'YLim'),'--','Linewidth',3,'color',matRad_cfg.gui.highlightColor);hold(handles.axesFig,'on');
                    
                else
                    PlotHandles{Cnt,1} =[];
                end
                
                Lines  = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),1);
                Labels = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),2);
                l=legend(handles.axesFig,[Lines{:}],Labels{:},'TextColor',matRad_cfg.gui.textColor,'FontSize',defaultFontSize,'Color',matRad_cfg.gui.backgroundColor,'EdgeColor',matRad_cfg.gui.textColor);
                xlabel('radiological depth [mm]','FontSize',matRad_cfg.gui.fontSize,'Color',matRad_cfg.gui.textColor);
                grid on, grid minor
                if ~isempty(rightAx)
                    set(rightAx,'XLim',get(handles.axesFig,'XLim'));
                end
            else
                % create legend for the visible VOI
                if this.typeOfPlot==2 || ~this.plotContour || isempty(this.AxesHandlesVOI) || isempty([this.AxesHandlesVOI{:}]) %isempty(find(this.VOIPlotFlag, 1))
                    l=legend(handles.axesFig,'off');
                else
%                    
                    % in case of multiple lines per VOI, only display the legend once
                    VOIlines=this.AxesHandlesVOI(~cellfun(@isempty,this.AxesHandlesVOI));
                    VOIlegendlines=cellfun(@(v)v(1),VOIlines);

                    l=legend(handles.axesFig, VOIlegendlines, this.cst{ this.sliceContourLegend,2}); %, 'FontSize',8
                    hold(handles.axesFig,'on');
                    
                end
            end
            set(l,'FontSize',defaultFontSize);
            if this.plotLegend
                set(l,'Visible','on')
            else
                set(l,'Visible','off')
            end
            this.legendHandle=l;
            
            
%             if this.rememberCurrAxes
%                 axis(handles.axesFig);%currAxes);
%             end
            
            hold(handles.axesFig,'off');
            
            %this.cBarChanged = false;
            
%             if this.typeOfPlot==1
%                 UpdateColormapOptions(handles);
%             end
            
            %this.legendHandle=legend(handles.axesFig);
            this.handles = handles;
            %profile off;
            %profile viewer;
            
            notifyPlotUpdated(this);
        end
       
        %Update IsodoseLines
        function this = updateIsoDoseLineCache(this)
            handles=this.handles;
            
            %Lock triggering an update during isoline caching
            currLock = this.updateLock;
            this.updateLock = true;
            
            if evalin('base','exist(''resultGUI'')')                             
                resultGUI = evalin('base','resultGUI');
                % select first cube if selected option does not exist
                if ~isfield(resultGUI,this.SelectedDisplayOption)
                    CubeNames = fieldnames(resultGUI);
                    cubeIx = structfun(@(s) numel(size(s)) == 3 && isnumeric(s),resultGUI);
                    CubeNames = CubeNames(cubeIx);
                    this.SelectedDisplayOption = CubeNames{1,1};
                else
                    
                    
                end
                dose = resultGUI.(this.SelectedDisplayOption);
                
                %if function is called for the first time then set display parameters                        
                if (isempty(this.dispWindow{2,2}) || ~this.checkUpdateNecessary({'viewer_options'},this.evt) ) &&  ~this.lockColorSettings 
                   this.dispWindow{2,1} = [min(dose(:)) max(dose(:))*1.001]; % set default dose range
                   this.dispWindow{2,2} = [min(dose(:)) max(dose(:))*1.001]; % set min max values
                end
                
                minMaxRange = this.dispWindow{2,1};
                % if upper colorrange is defined then use it otherwise 120% iso dose
                upperMargin = 1;
                if abs((max(dose(:)) - this.dispWindow{2,1}(1,2))) < 0.01  * max(dose(:))
                    upperMargin = 1.2;
                end
                
                %this creates a loop(needed the first time a dose cube is loaded)
                if isempty(this.IsoDose_Levels) || ~this.NewIsoDoseFlag || ~this.checkUpdateNecessary({'viewer_options'},this.evt)
                    vLevels                  = [0.1:0.1:0.9 0.95:0.05:upperMargin];
                    referenceDose            = (minMaxRange(1,2))/(upperMargin);
                    
                    this.IsoDose_Levels   = minMaxRange(1,1) + (referenceDose-minMaxRange(1,1)) * vLevels;

                end
                
                
                this.IsoDose_Contours = matRad_computeIsoDoseContours(dose,this.IsoDose_Levels);
            end
            this.handles = handles;
            
            this.updateLock = currLock;
        end
        
        %% Data Cursors
        function cursorText = dataCursorUpdateFunction(this,hObject, eventdata)%event_obj)
            % Display the position of the data cursor
            % obj          Currently not used (empty)
            % event_obj    Handle to event object
            % output_txt   Data cursor text string (string or cell array of strings).
            
            target = findall(0,'Name','matRadGUI');
            
            % Get GUI data (maybe there is another way?)
            %handles = guidata(target);
            
            % position of the data point to label
            pos = get(eventdata,'Position');
            
            %Different behavior for image and profile plot
            if this.typeOfPlot ==1 %get(handles.popupTypeOfPlot,'Value')==1 %Image view
                cursorText = cell(0,1);
                try
                    
                    if evalin('base','exist(''ct'')') %handles.State >= 1
%                         plane = get(handles.popupPlane,'Value');
%                         slice = round(get(handles.sliderSlice,'Value'));
                        
                        %Get the CT values
                        ct  = evalin('base','ct');
                        
                        %We differentiate between pos and ix, since the user may put
                        %the datatip on an isoline which returns a continous position
                        cubePos = zeros(1,3);
                        cubePos(this.plane) = this.slice;
                        cubePos(1:end ~= this.plane) = fliplr(pos);
                        cubeIx = round(cubePos);
                        %Here comes the index permutation stuff
                        %Cube Index
                        cursorText{end+1,1} = ['Cube Index: ' mat2str(cubeIx)];
                        %Space Coordinates
                        coords = matRad_cubeIndex2worldCoords(cubeIx,ct);
                        cursorText{end+1,1} = ['Space Coordinates: ' mat2str(coords,5) ' mm'];
                        
                        ctVal = ct.cubeHU{1}(cubeIx(1),cubeIx(2),cubeIx(3));
                        cursorText{end+1,1} = ['HU Value: ' num2str(ctVal,3)];
                    end
                catch
                    cursorText{end+1,1} = 'Error while retreiving CT Data!';
                end
                
                
                    %Add dose information if available
                    if evalin('base','exist(''resultGUI'')') %handles.State == 3
                        %get result structure
                        result = evalin('base','resultGUI');
                        
                        %get all cubes from the ResultGUI 
                        resultNames = fieldnames(result); %get(handles.popupDisplayOption,'String');
                        
                        %Display all values of fields found in the resultGUI struct
                        for runResult = 1:numel(resultNames)
                            if ~isstruct(result.(resultNames{runResult,1})) && ~isvector(result.(resultNames{runResult,1}))
                            %try
                                name = resultNames{runResult};
                                if isfield(result,name) % (check the dimensions, same as CT)
                                    field = result.(name);
                                    val = field(cubeIx(1),cubeIx(2),cubeIx(3));
                                    cursorText{end+1,1} = [name ': ' num2str(val,3)];
                                end
                            %catch
                                %cursorText{end+1,1} = 'Error while retreiving Data!';
                            end
                        end
                    end
                
            else %Profile view
                cursorText = cell(2,1);
                cursorText{1} = ['Radiological Depth: ' num2str(pos(1),3) ' mm'];
                cursorText{2} = [get(target,'DisplayName') ': ' num2str(pos(2),3)];
            end
            
        end
        
        %Scroll wheel update
        function matRadScrollWheelFcn(this,src,event)
            % Check Position
            cursorPos = get(src,'CurrentPoint');
            viewerPos = get(this.widgetHandle,'Position');

            %If the request came from within the widget change the slice
            if cursorPos(1) > viewerPos(1) && cursorPos(1) < viewerPos(1) + viewerPos(3) ...
                    && cursorPos(2) > viewerPos(2) && cursorPos(2) < viewerPos(2) + viewerPos(4)

                % compute new slice
                this.slice= this.slice - event.VerticalScrollCount;
            end
        end
        
        %Toggle Legend
        function legendToggleFunction(this,src,event)
            if isempty(this.legendHandle) || ~isobject(this.legendHandle)
                return;
            end
            if this.plotLegend
                set(this.legendHandle,'Visible','on')
            else
                set(this.legendHandle,'Visible','off')
            end
        end
        
        %Toggle Colorbar 
        function colorBarToggleFunction(this,src,event)
            if isempty(this.cBarHandle) || ~isobject(this.cBarHandle) || this.updateLock
                return;
            end
            if this.plotColorBar
                
                if evalin('base','exist(''resultGUI'')')
                    this.colorData=2;
                else evalin('base','exist(''ct'')')
                    this.colorData=1;
                end
                set(this.cBarHandle,'Visible','on')
            else
                set(this.cBarHandle,'Visible','off');
            end            
            % send a notification that the plot has changed (to update the options)
            %this.notifyPlotUpdated();
        end
        
        %
        function initValues(this)
            lockState=this.updateLock;
            
            if lockState 
                return;
            end
            
            this.updateLock=true;
            
            if isempty(this.plane)
                this.plane=3;
            end
            
            if evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')                 
                % update slice, beam and offset sliders parameters
                
                ct = evalin('base','ct');
                cst = evalin('base','cst');
                cst = matRad_computeVoiContoursWrapper(cst,ct);
                assignin('base','cst',cst);
                this.cst = cst;
                
                % define context menu for structures
                this.VOIPlotFlag=false(size(this.cst,1),1);
                for i = 1:size(this.cst,1)
                    if this.cst{i,5}.Visible
                        this.VOIPlotFlag(i) = true;
                    end
                end

                if isfield(ct, 'cubeHU')
                    minMax = [min(ct.cubeHU{1}(:)) max(ct.cubeHU{1}(:))];

                    if diff(minMax) == 0
                        minMax = [-1000 2000];
                    end
                else
                    minMax = [min(ct.cube{1}(:)) max(ct.cube{1}(:))];

                    if diff(minMax) == 0
                        minMax = [0 2];
                    end
                end


                planeCenters = ceil(ct.cubeDim./ 2);
                this.numOfBeams = 1;

                visQuantity = this.tryVisQuantityFromPln();
                
                if evalin('base','exist(''pln'')')
                    pln = evalin('base','pln');
                    if isfield(pln,'propStf') && isfield(pln.propStf,'isoCenter')
                        isoCoordinates = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:), ct);
                        planeCenters = ceil(isoCoordinates);
                        this.numOfBeams=pln.propStf.numOfBeams;
                    end
                end
                       
                this.slice = planeCenters(this.plane);                        

                % set profile offset slider
                this.OffsetMinMax = [-100 100];
                vRange = sum(abs(this.OffsetMinMax));

                if strcmp(this.ProfileType,'lateral')
                    this.OffsetSliderStep = vRange/ct.resolution.x;
                else
                    this.OffsetSliderStep = vRange/ct.resolution.y;
                end
                this.OffsetSliderStep=[1/this.OffsetSliderStep 1/this.OffsetSliderStep];

                               
                selectionIndex=1;
                this.plotColorBar=true;
               
                
                if evalin('base','exist(''resultGUI'')')
                    this.colorData=2;
                    this.plotColorBar=true;
                    selectionIndex=1;
                    
                    Result = evalin('base','resultGUI');
                    
                    this.DispInfo = fieldnames(Result);
                    
                    this.updateDisplaySelection(visQuantity);
                else
                    this.colorData=1;
                    this.SelectedDisplayAllOptions = 'no option available';
                    this.SelectedDisplayOption = '';
                end
            else %no data is loaded 
                this.slice=1;
                this.numOfBeams=1;
                this.OffsetMinMax = [1 1];
                this.profileOffset=1;
                this.OffsetSliderStep=[1 1];
                this.colorData=1;
                this.plotColorBar=false;
                selectionIndex=1;
                minMax = [0 1];
                this.SelectedDisplayAllOptions = 'no option available';
                this.SelectedDisplayOption = '';
            end
            
            this.dispWindow{selectionIndex,1} = minMax;
            this.dispWindow{selectionIndex,2} = minMax;
             
            this.updateLock=lockState;
        end
        
        %update the Viewer
        function updateValues(this)
            lockState=this.updateLock;
            
            if lockState
                return;
            end
            
            this.updateLock=true;
                                   
            if evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')
                % update slice, beam and offset sliders parameters
                ct = evalin('base','ct');
                cst = evalin('base','cst');
                this.cst = cst;
                
                % define context menu for structures
                this.VOIPlotFlag=false(size(this.cst,1),1);
                for i = 1:size(this.cst,1)
                    if this.cst{i,5}.Visible
                        this.VOIPlotFlag(i) = true;
                    end
                end
                % set isoCenter values 
                % Note: only defined for the first Isocenter
                if evalin('base','exist(''pln'')')
                    pln = evalin('base','pln');
                    if isfield(pln,'propStf') && isfield(pln.propStf,'isoCenter')
                        this.vIsoCenter      = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:), ct);
                    else
                        this.plotIsoCenter = false;
                    end
                else
                    this.plotIsoCenter = false;
                end

                 % set profile offset slider
                this.OffsetMinMax = [-100 100];
                vRange = sum(abs(this.OffsetMinMax));
                
                if strcmp(this.ProfileType,'lateral')
                    this.OffsetSliderStep = vRange/ct.resolution.x;
                else
                    this.OffsetSliderStep = vRange/ct.resolution.y;
                end
                this.OffsetSliderStep=[1/this.OffsetSliderStep 1/this.OffsetSliderStep];

                this.updateDisplaySelection();
            end

            this.updateLock=lockState;
        end
    
        function updateDisplaySelection(this,visSelection)
            %Lock triggering an update during isoline caching
            currLock = this.updateLock;
            this.updateLock = true;

            if nargin < 2
                visSelection = [];
            end

            if evalin('base','exist(''resultGUI'')')                   
                result = evalin('base','resultGUI');
                
                this.DispInfo = fieldnames(result);
                for i = 1:size(this.DispInfo,1)
                    
                    % delete weight vectors in Result struct for plotting
                    if isstruct(result.(this.DispInfo{i,1})) || isvector(result.(this.DispInfo{i,1}))
                        result = rmfield(result,this.DispInfo{i,1});
                        this.DispInfo{i,2}=false;
                    else
                        %second dimension indicates if it should be plotted
                        this.DispInfo{i,2} = true;
                        % determine units (third dimension) and left or
                        % right axis (fourth dimension, ignored so far
                        if strfind(this.DispInfo{i,1},'physicalDose')
                            this.DispInfo{i,3} = 'Gy';
                            this.DispInfo{i,4} = 'left';   
                        elseif strfind(this.DispInfo{i,1},'alpha')
                            this.DispInfo{i,3} = 'Gy^{-1}';
                            this.DispInfo{i,4} = 'left';   
                        elseif strfind(this.DispInfo{i,1},'beta')
                            this.DispInfo{i,3} = 'Gy^{-2}';
                            this.DispInfo{i,4} = 'left';   
                        elseif strfind(this.DispInfo{i,1},'RBExDose')
                            this.DispInfo{i,3} = 'Gy(RBE)';
                            this.DispInfo{i,4} = 'left';   
                        elseif strfind(this.DispInfo{i,1},'LET')
                            this.DispInfo{i,3} = 'keV/um';
                            this.DispInfo{i,4} = 'left';   
                        elseif strfind(this.DispInfo{i,1},'effect')
                            this.DispInfo{i,3} = '1';
                            this.DispInfo{i,4} = 'right';   
                        elseif strfind(this.DispInfo{i,1},'RBE')
                            this.DispInfo{i,3} = '1';
                            this.DispInfo{i,4} = 'right';   
                        elseif strfind(this.DispInfo{i,1},'BED')
                            this.DispInfo{i,3} = 'Gy';
                            this.DispInfo{i,4} = 'left';   
                        else
                            this.DispInfo{i,3} = 'a.u.';
                            this.DispInfo{i,4} = 'right';
                        end
                    end
                end
                
                this.SelectedDisplayAllOptions=fieldnames(result);                    
                
                if ~isempty(visSelection) && isfield(result,visSelection)
                    this.SelectedDisplayOption = visSelection;
                elseif ~isfield(result,this.SelectedDisplayOption)
                    this.SelectedDisplayOption = this.tryVisQuantityFromPln('physicalDose');
                else
                    %Keep option
                end
                
                if ~any(strcmp(this.SelectedDisplayOption,fieldnames(result)))
                    this.SelectedDisplayOption = this.tryVisQuantityFromPln('physicalDose');
                    if ~any(strcmp(this.SelectedDisplayOption,fieldnames(result)))
                        this.SelectedDisplayOption = this.DispInfo{find([this.DispInfo{:,2}],1,'first'),1};
                    end

                    this.updateIsoDoseLineCache();
                end               
            else
                this.SelectedDisplayAllOptions = 'no option available';
                this.SelectedDisplayOption = '';
            end

            this.updateLock = currLock;
        end        

        function visQuantity = tryVisQuantityFromPln(~, default)
            if nargin < 2
                default = [];
            end
            visQuantity = default;
            if evalin('base','exist(''pln'')')
                pln = evalin('base','pln');
                if isfield(pln,'propOpt') && isfield(pln.propOpt, 'quantityOpt')
                    switch pln.propOpt.quantityOpt
                        case 'physicalDose'
                            visQuantity = 'physicalDose';
                        case {'RBExDose', 'effect'}
                            visQuantity = 'RBExDose';
                        otherwise
                            %Do Nothing
                    end
                end
            end
        end
    end

    
end