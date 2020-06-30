classdef matRad_VisualizationWidget < matRad_Widget
    
    properties
        viewingWidgetHandle;
    end
    
    methods
        function this = matRad_VisualizationWidget(handleParent,viewingWidgetHandle)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[170 45 80 10],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],... 'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Visualization',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
                
            end
            this = this@matRad_Widget(handleParent);  
            
            handles=this.handles;
            
            if nargin==2
                this.viewingWidgetHandle=viewingWidgetHandle;
            else
                set(handles.btnDVH,'Enable','off');
                set(handles.popupDisplayOption,'Enable','off');
                set(handles.btnProfileType,'Enable','off');
                set(handles.popupTypeOfPlot,'Enable','off');
                set(handles.popupPlane,'Enable','off');
                set(handles.radiobtnCT,'Enable','off');
                set(handles.radiobtnContour,'Enable','off');
                set(handles.radiobtnDose,'Enable','off');
                set(handles.radiobtnIsoDoseLines,'Enable','off');
                set(handles.sliderSlice,'Enable','off');
                set(handles.radiobtnIsoDoseLinesLabels,'Enable','off');
                set(handles.radioBtnIsoCenter,'Enable','off');
                set(handles.radiobtnPlan,'Enable','off');
                set(handles.btn3Dview,'Enable','off');
            end
            this.handles=handles;
        end
        
        
        function this = initialize(this)
            
        end
        
        
        
        function this = update(this)
            if isa(this.viewingWidgetHandle,'matRad_ViewingWidget')                
                % get the default values from the viewer widget
                this.getFromViewingWidget();
            else
                handles=this.handles;
                % disable all buttons
                set(handles.popupDisplayOption,'Enable','off');
                set(handles.btnProfileType,'Enable','off');
                set(handles.popupTypeOfPlot,'Enable','off');
                set(handles.popupPlane,'Enable','off');
                set(handles.radiobtnCT,'Enable','off');
                set(handles.radiobtnContour,'Enable','off');
                set(handles.radiobtnDose,'Enable','off');
                set(handles.radiobtnIsoDoseLines,'Enable','off');
                set(handles.sliderSlice,'Enable','off');
                set(handles.radiobtnIsoDoseLinesLabels,'Enable','off');
                set(handles.radioBtnIsoCenter,'Enable','off');
                set(handles.radiobtnPlan,'Enable','off');
                this.handles=handles;
            end
        end
                
%         function viewingWidgetHandle=get.viewingWidgetHandle(this)
%             viewingWidgetHandle=this.viewingWidgetHandle;
%         end
        
        function set.viewingWidgetHandle(this,value)
            if isa(value,'matRad_ViewingWidget')
                this.viewingWidgetHandle=value;
                
                % get the default values from the viewer widget      
                this.getFromViewingWidget();
                
            else
                handles=this.handles;
                % disable all buttons
                set(handles.popupDisplayOption,'Enable','off');
                set(handles.btnProfileType,'Enable','off');
                set(handles.popupTypeOfPlot,'Enable','off');
                set(handles.popupPlane,'Enable','off');
                set(handles.radiobtnCT,'Enable','off');
                set(handles.radiobtnContour,'Enable','off');
                set(handles.radiobtnDose,'Enable','off');
                set(handles.radiobtnIsoDoseLines,'Enable','off');
                set(handles.sliderSlice,'Enable','off');
                set(handles.radiobtnIsoDoseLinesLabels,'Enable','off');
                set(handles.radioBtnIsoCenter,'Enable','off');
                set(handles.radiobtnPlan,'Enable','off');
                this.handles=handles;
            end
        end
    end
    
    methods (Access = protected)
        function this = createLayout(this)
            h36 = this.widgetHandle;
            h37 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String',{  'coronal'; 'sagital'; 'axial' },...
                'Style','popupmenu',...
                'Value',3,...
                'Position',[0.465315808689303 0.582191780821918 0.113636363636364 0.143835616438356],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) popupPlane_Callback(this,hObject,eventdata),...
                'Tag','popupPlane',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h38 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String',{  'Slider' },...
                'Style','slider',...
                'Position',[0.134961439588689 0.796610169491525 0.167095115681234 0.096045197740113],...
                'BackgroundColor',[0.9 0.9 0.9],...
                'Callback',@(hObject,eventdata) sliderSlice_Callback(this,hObject,eventdata),...
                'BusyAction','cancel',...
                'Interruptible','off',...
                'FontSize',8,...
                'Tag','sliderSlice');
            
            h39 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','Plane Selection',...
                'Style','text',...
                'Position',[0.32 0.55 0.11969696969697 0.2],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtPlanSelection');
            
            h40 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Slice Selection',...
                'Style','text',...
                'Position',[0.00909090909090909 0.75 0.113636363636364 0.2],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','text9' );
            
            h41 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot contour',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.780445969125214 0.733212341197822 0.169811320754717 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnContour_Callback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','radiobtnContour' );
            
            h42 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot dose',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.780445969125214 0.466969147005444 0.169811320754717 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnDose_Callback(this,hObject,eventdata),...
                'Children',[],...
                'FontSize',8,...
                'Tag','radiobtnDose' );
            
            
            h43 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot isolines',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.780445969125214 0.600907441016334 0.150943396226415 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnIsoDoseLines_Callback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','radiobtnIsoDoseLines');
            
            h44 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Type of plot',...
                'Style','text',...
                'Position',[0.32 0.793103448275862 0.2 0.124137931034483],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtTypeOfPlot' );
            
            h45 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String',{  'intensity'; 'profile' },...
                'Style','popupmenu',...
                'Value',1,...
                'Position',[0.465315808689303 0.801369863013699 0.113636363636364 0.143835616438356],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata) popupTypeOfPlot_Callback(this,hObject,eventdata),...
                'Tag','popupTypeOfPlot',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h46 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Display option',...
                'Style','text',...
                'Position',[0.32 0.36 0.2 0.102739726027397],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtDisplayOption' );
            
            h47 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','Please select ... ',...
                'Style','popupmenu',...
                'Value',1,...
                'Position',[0.465315808689303 0.36986301369863 0.196969696969697 0.136986301369863],...
                'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
                'Callback',@(hObject,eventdata)popupDisplayOption_Callback(this,hObject,eventdata),...
                'Tag','popupDisplayOption',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h48 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Beam Selection',...
                'Style','text',...
                'Position',[0.00857632933104631 0.5 0.118353344768439 0.2],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','txtBeamSelection' );
            
            h49 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','SliderBeamSelection',...
                'Style','slider',...
                'Position',[0.134961439588689 0.542372881355932 0.167095115681234 0.096045197740113],...
                'BackgroundColor',[0.9 0.9 0.9],...
                'Callback',@(hObject,eventdata) sliderBeamSelection_Callback(this,hObject,eventdata),...
                'Enable','off',...
                'FontSize',8,...
                'Tag','sliderBeamSelection');
            
            h50 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','lateral',...
                'Position',[0.658025928757913 0.794520547945205 0.0863636363636364 0.157534246575343],...
                'BackgroundColor',[0.8 0.8 0.8],...
                'Callback',@(hObject,eventdata) btnProfileType_Callback(this,hObject,eventdata),...
                'Enable','off',...
                'Tag','btnProfileType',...
                'FontSize',8,...
                'FontWeight','bold' );
            
            h51 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','GoTo',...
                'Style','text',...
                'Position',[0.59 0.801369863013698 0.06 0.123287671232877],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Children',[],...
                'FontSize',8,...
                'Tag','text16' );
            
            h52 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','Show DVH/QI',...
                'Position',[0.51413881748072 0.0677966101694915 0.153393316195373 0.129943502824859],...
                'BackgroundColor',[0.8 0.8 0.8],...
                'Callback',@(hObject,eventdata) btnDVH_Callback(this,hObject,eventdata),...
                'Enable','off',...
                'Tag','btnDVH',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h53 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot isolines labels',...
                'Style','radiobutton',...
                'Position',[0.780445969125214 0.343557168784029 0.252401372212693 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnIsoDoseLinesLabels_Callback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','radiobtnIsoDoseLinesLabels');
            
            h54 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Offset',...
                'Style','text',...
                'Position',[0.00909090909090909 0.25 0.118181818181818 0.123287671232877],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'FontSize',8,...
                'Tag','textOffset');
            
            h55 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','SliderOffset',...
                'Style','slider',...
                'Position',[0.134961439588689 0.271186440677966 0.167095115681234 0.096045197740113],...
                'BackgroundColor',[0.9 0.9 0.9],...
                'Callback',@(hObject,eventdata) sliderOffset_Callback(this,hObject,eventdata),...
                'Enable','off',...
                'FontSize',8,...
                'Tag','sliderOffset');
            
            
            h56 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot iso center',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.780445969125214 0.205989110707804 0.2 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radioBtnIsoCenter_Callback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','radioBtnIsoCenter');
            
            
            h57 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','Open 3D-View',...
                'Position',[0.595848595848596 0.578947368421053 0.148962148962149 0.157894736842105],...
                'BackgroundColor',[0.8 0.8 0.8],...
                'Callback',@(hObject,eventdata) btn3Dview_Callback(this,hObject,eventdata),...
                'Enable','off',...
                'Tag','btn3Dview',...
                'FontSize',8,...
                'FontWeight','bold');
            
            h58 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot CT',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.780445969125214 0.864791288566243 0.169811320754717 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnCT_Callback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','radiobtnCT');
            
            
            h59 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','visualize plan/beams',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.78021978021978 0.0736842105263158 0.3002442002442 0.115789473684211],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnPlan_Callback(this,hObject,eventdata),...
                'FontSize',8,...
                'Tag','radiobtnPlan');
            
            
            this.createHandles();
        end
    end
    
    methods (Access = protected)
        function getFromViewingWidget(this)
            handles=this.handles;
            if strcmp(this.viewingWidgetHandle.ProfileType,'lateral')
                set(handles.btnProfileType,'String','longitudinal');
            else
                set(handles.btnProfileType,'String','lateral');
            end
            
            set(handles.popupTypeOfPlot,'Value',this.viewingWidgetHandle.typeOfPlot);
            set(handles.popupPlane,'Value',this.viewingWidgetHandle.plane);
            set(handles.radiobtnContour,'Value',this.viewingWidgetHandle.plotContour);
            set(handles.radiobtnDose,'Value',this.viewingWidgetHandle.plotDose);
            set(handles.radiobtnIsoDoseLines,'Value',this.viewingWidgetHandle.plotIsoDoseLines);
            set(handles.radiobtnIsoDoseLinesLabels,'Value',this.viewingWidgetHandle.plotIsoDoseLinesLabels);
            set(handles.radioBtnIsoCenter,'Value',this.viewingWidgetHandle.plotIsoCenter);
            set(handles.radiobtnPlan,'Value',this.viewingWidgetHandle.plotPlan);
            
            % update the sliders 
            set(handles.sliderSlice,'Min',1,'Max',this.viewingWidgetHandle.maxSlice,...
                'Value', this.viewingWidgetHandle.slice, ...
                'SliderStep',this.viewingWidgetHandle.SliceSliderStep);
            
            if this.viewingWidgetHandle.numofBeams>1
                set(handles.sliderBeamSelection,'Min',1,'Max',this.viewingWidgetHandle.numOfBeams,...
                    'Value',this.viewingWidgetHandle.selectedBeam,...
                    'SliderStep',[1/(this.viewingWidgetHandle.numOfBeams-1) 1/(this.viewingWidgetHandle.numOfBeams-1)]);
            else
                set(handles.sliderBeamSelection,'Min',1,'Max',1, 'Value',1,'SliderStep',[1 1]);
            end
           
            set(handles.sliderOffset,'Min',this.viewingWidgetHandle.OffsetMinMax(1),'Max',this.viewingWidgetHandle.OffsetMinMax(2),...
                'Value',this.viewingWidgetHandle.profileOffset,...
                'SliderStep',this.viewingWidgetHandle.OffsetSliderStep);
            
            set(handles.popupDisplayOption,'String',this.viewingWidgetHandle.SelectedDisplayAllOptions);
            if ~strcmp(this.viewingWidgetHandle.SelectedDisplayOption,'')
                set(handles.popupDisplayOption,'Value',find(strcmp(this.viewingWidgetHandle.SelectedDisplayOption,this.viewingWidgetHandle.SelectedDisplayAllOptions)));
            end
            
            if strcmp(this.viewingWidgetHandle.SelectedDisplayOption,'') % no data is loaded
                % disable 3D and DVH button
                set(handles.btn3Dview,'Enable','off');
                set(handles.btnDVH,'Enable','off');
            else
                set(handles.btn3Dview,'Enable','on');
                
                if evalin('base','exist(''resultGUI'')')
                    set(handles.btnDVH,'Enable','on');
                else
                    set(handles.btnDVH,'Enable','off');
                end
                
                 %% enable and diasble buttons according to type of plot
                % intensity plot
                if this.viewingWidgetHandle.typeOfPlot == 1
                    
                    set(handles.sliderBeamSelection,'Enable','off')
                    set(handles.sliderOffset,'Enable','off')
                    set(handles.popupDisplayOption,'Enable','on')
                    set(handles.btnProfileType,'Enable','off');
                    set(handles.popupPlane,'Enable','on');
                    set(handles.radiobtnCT,'Enable','on');
                    set(handles.radiobtnContour,'Enable','on');
                    set(handles.radiobtnDose,'Enable','on');
                    set(handles.radiobtnIsoDoseLines,'Enable','on');
                    set(handles.radiobtnIsoDoseLinesLabels,'Enable','on');
                    set(handles.sliderSlice,'Enable','on');
                    
                    % profile plot
                elseif this.viewingWidgetHandle.typeOfPlot == 2
                    
                    set(handles.popupDisplayOption,'Enable','on');
                    set(handles.btnProfileType,'Enable','on');
                    set(handles.popupPlane,'Enable','off');
                    set(handles.radiobtnCT,'Enable','off');
                    set(handles.radiobtnContour,'Enable','off');
                    set(handles.radiobtnDose,'Enable','off');
                    set(handles.radiobtnIsoDoseLines,'Enable','off');
                    set(handles.sliderSlice,'Enable','off');
                    set(handles.radiobtnIsoDoseLinesLabels,'Enable','off');
                    set(handles.btnProfileType,'Enable','on')
                    
                end
                
            end
            this.handles=handles;
        end        
        
         % H37 Calback
        function popupPlane_Callback(this, hObject, event)
            % hObject    handle to popupPlane (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: contents = cellstr(get(hObject,'String')) returns popupPlane contents as cell array
            %        contents{get(hObject,'Value')} returns selected item from popupPlane
            
            % set slice slider
            handles = this.handles;
            
            this.viewingWidgetHandle.plane = get(handles.popupPlane,'value');
%             try
%                 if evalin('base','exist(''pln'')') && evalin('base','exist(''ct'')') && ...
%                         evalin('base','exist(''cst'')') %if handles.State > 0
%                     ct = evalin('base', 'ct');
%                     set(handles.sliderSlice,'Min',1,'Max',ct.cubeDim(handles.plane),...
%                         'SliderStep',[1/(ct.cubeDim(handles.plane)-1) 1/(ct.cubeDim(handles.plane)-1)]);
%                     if ~evalin('base','exist(''ResultGUI'')') %handles.State < 3
%                         set(handles.sliderSlice,'Value',round(ct.cubeDim(handles.plane)/2));
%                     else
%                         pln = evalin('base','pln');
%                         
%                         if handles.plane == 1
%                             set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,2)/ct.resolution.x));
%                         elseif handles.plane == 2
%                             set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,1)/ct.resolution.y));
%                         elseif handles.plane == 3
%                             set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,3)/ct.resolution.z));
%                         end
%                         
%                     end
%                 end
%             catch
%             end
             
%             handles.rememberCurrAxes = false;
%             UpdatePlot(handles);
%             handles.rememberCurrAxes = true;
            this.handles = handles;
        end
        
         %45 Callback
        function popupTypeOfPlot_Callback(this, hObject, event)
            this.viewingWidgetHandle.typeOfPlot=get(hObject,'Value');
            handles = this.handles;
            
            % intensity plot
            if get(hObject,'Value') == 1
                
                set(handles.sliderBeamSelection,'Enable','off')
                set(handles.sliderOffset,'Enable','off')
                set(handles.popupDisplayOption,'Enable','on')
                set(handles.btnProfileType,'Enable','off');
                set(handles.popupPlane,'Enable','on');
                set(handles.radiobtnCT,'Enable','on');
                set(handles.radiobtnContour,'Enable','on');
                set(handles.radiobtnDose,'Enable','on');
                set(handles.radiobtnIsoDoseLines,'Enable','on');
                set(handles.radiobtnIsoDoseLinesLabels,'Enable','on');
                set(handles.sliderSlice,'Enable','on');
                
                % profile plot
            elseif get(hObject,'Value') == 2   
                     
                if evalin('base','exist(''pln'')') && evalin('base','exist(''ct'')')
                    if this.viewingWidgetHandle.numofBeams>1
                        set(handles.sliderBeamSelection,'Enable','on');
                    end
                    set(handles.sliderOffset,'Enable','on');
                end
                                
                set(handles.popupDisplayOption,'Enable','on');
                set(handles.btnProfileType,'Enable','on');
                set(handles.popupPlane,'Enable','off');
                set(handles.radiobtnCT,'Enable','off');
                set(handles.radiobtnContour,'Enable','off');
                set(handles.radiobtnDose,'Enable','off');
                set(handles.radiobtnIsoDoseLines,'Enable','off');
                set(handles.sliderSlice,'Enable','off');
                set(handles.radiobtnIsoDoseLinesLabels,'Enable','off');
                  
                set(handles.btnProfileType,'Enable','on')
%                 
%                 if strcmp(get(handles.btnProfileType,'String'),'lateral')
%                     this.viewingWidgetHandle.ProfileType = 'longitudinal';
%                 else
%                     this.viewingWidgetHandle.ProfileType = 'lateral';
%                 end
            end
            
            %this.viewingWidgetHandle.cBarChanged = true;
            
%             handles.rememberCurrAxes = false;
%             cla(handles.axesFig,'reset');
%             UpdatePlot(handles);
%             handles.rememberCurrAxes = true;
            
            this.handles = handles;
        end
        
          % 47 Callback
        function popupDisplayOption_Callback(this, hObject, event)
            %this.updateIsodoseLine();
            content = get(hObject,'String');
            if strcmp(content,'no option available')
                return
            end
            
            handles = this.handles;
            this.viewingWidgetHandle.SelectedDisplayOption = content{get(hObject,'Value'),1};
            %this.viewingWidgetHandle.SelectedDisplayOptionIdx = get(hObject,'Value');
            
%             if ~isfield(handles,'colormapLocked') || ~handles.colormapLocked
%                 handles.dispWindow{3,1} = []; handles.dispWindow{3,2} = [];
%             end
            
            %this.viewingWidgetHandle.updateIsoDoseLineCache();
            
%             if evalin('base','exist(''resultGUI'')')
%                 
%                 resultGUI = evalin('base','resultGUI');
%                 % select first cube if selected option does not exist
%                 if ~isfield(resultGUI,this.viewingWidgetHandle.SelectedDisplayOption)
%                     CubeNames = fieldnames(resultGUI);
%                     dose = resultGUI.(CubeNames{1,1});
%                 else
%                     dose = resultGUI.(this.viewingWidgetHandle.SelectedDisplayOption);
%                 end
%                 
%                 %if function is called for the first time then set display parameters
%                 if isempty(this.viewingWidgetHandle.dispWindow{3,2})
%                     this.viewingWidgetHandle.dispWindow{3,1} = [min(dose(:)) max(dose(:))]; % set default dose range
%                     this.viewingWidgetHandle.dispWindow{3,2} = [min(dose(:)) max(dose(:))]; % set min max values
%                 end
%                 
%                 minMaxRange = this.viewingWidgetHandle.dispWindow{3,1};
%                 % if upper colorrange is defined then use it otherwise 120% iso dose
%                 upperMargin = 1;
%                 if abs((max(dose(:)) - this.viewingWidgetHandle.dispWindow{3,1}(1,2))) < 0.01  * max(dose(:))
%                     upperMargin = 1.2;
%                 end
%                 
%                 %if (length(this.viewingWidgetHandle.IsoDose_Levels) == 1 && this.viewingWidgetHandle.IsoDose_Levels(1,1) == 0) || ~this.NewIsoDoseFlag
%                     vLevels                  = [0.1:0.1:0.9 0.95:0.05:upperMargin];
%                     referenceDose            = (minMaxRange(1,2))/(upperMargin);
%                     this.viewingWidgetHandle.IsoDose_Levels   = minMaxRange(1,1) + (referenceDose-minMaxRange(1,1)) * vLevels;
%                 %end 
%             end
            
            %this.viewingWidgetHandle.cBarChanged = true;
            
            %UpdatePlot(handles);
            this.handles = handles;
            
            
            
        end
        
        % H49 Callback
        function sliderBeamSelection_Callback(this, hObject, event)
            % hObject    handle to sliderBeamSelection (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            
            handles = this.handles;
            this.viewingWidgetHandle.selectedBeam = round(get(hObject,'Value'));
            set(hObject, 'Value', this.viewingWidgetHandle.selectedBeam);
%             handles.rememberCurrAxes = false;
%             UpdatePlot(handles);
%             handles.rememberCurrAxes = true;
            
            this.handles = handles;
        end
         
         % 50 Callback
        function btnProfileType_Callback(this, hObject, event)
            % hObject    handle to btnProfileType (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            %handles = this.handles;
            
            if strcmp(get(hObject,'Enable') ,'on')
                if strcmp(this.viewingWidgetHandle.ProfileType,'lateral')
                    this.viewingWidgetHandle.ProfileType = 'longitudinal';
                    set(hObject,'String','lateral');
                else
                    this.viewingWidgetHandle.ProfileType = 'lateral';
                    set(hObject,'String','longitudinal');
                end
                
%                 handles.rememberCurrAxes = false;
%                 UpdatePlot(handles);
%                 handles.rememberCurrAxes = true;
            end
            %this.handles = handles;
        end
         
        % 52 Callback
        function btnDVH_Callback(this, hObject, event)
            matRad_DVHWidget;
%             handles = this.handles;
%             resultGUI = evalin('base','resultGUI');
%             Content = get(handles.popupDisplayOption,'String');
%             
%             SelectedCube = Content{get(handles.popupDisplayOption,'Value')};
%             
%             pln = evalin('base','pln');
%             resultGUI_SelectedCube.physicalDose = resultGUI.(SelectedCube);
%             
%             if ~strcmp(pln.propOpt.bioOptimization,'none')
%                 
%                 %check if one of the default fields is selected
%                 if sum(strcmp(SelectedCube,{'physicalDose','effect','RBE,','RBExDose','alpha','beta'})) > 0
%                     resultGUI_SelectedCube.physicalDose = resultGUI.physicalDose;
%                     resultGUI_SelectedCube.RBExDose     = resultGUI.RBExDose;
%                 else
%                     Idx    = find(SelectedCube == '_');
%                     SelectedSuffix = SelectedCube(Idx(1):end);
%                     resultGUI_SelectedCube.physicalDose = resultGUI.(['physicalDose' SelectedSuffix]);
%                     resultGUI_SelectedCube.RBExDose     = resultGUI.(['RBExDose' SelectedSuffix]);
%                 end
%             end
%             
%             %adapt visibilty
%             cst = evalin('base','cst');
% %             for i = 1:size(cst,1)
% %                 cst{i,5}.Visible = this.viewingWidgetHandle.VOIPlotFlag(i);
% %             end
%             DVHfig=figure;
%             matRad_indicatorWrapper(DVHfig,cst,pln,resultGUI_SelectedCube);
%             
% %             assignin('base','cst',cst);
%             this.handles = handles;
            
        end
        
        %H55 Callback
        function sliderOffset_Callback(this, hObject, event)
            this.viewingWidgetHandle.profileOffset = get(hObject,'Value');
            %UpdatePlot(handles);
        end
       
        % 57 Callback
        function btn3Dview_Callback(this,hObject, event)
            
            matRad_3DWidget(this.viewingWidgetHandle);
            
            % hObject    handle to btn3Dview (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
%             handles = this.handles;
%             
%             if ~isfield(handles,'axesFig3D') || ~isfield(handles,'axesFig3D') || ~isgraphics(handles.axesFig3D)
%                 handles.fig3D = figure('Name','matRad 3D View');
%                 handles.axesFig3D = axes('Parent',handles.fig3D);
%                 view(handles.axesFig3D,3);
%             end
%             
%             %UpdatePlot(handles);
%             this.handles = handles;
%             this.viewingWidgetHandle.UpdatePlot();
        end
       % --- Executes on button press in radiobtnContour.
       function radiobtnContour_Callback(this,hObject, ~)
           % hObject    handle to radiobtnContour (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           
           % Hint: get(hObject,'Value') returns toggle state of radiobtnContour
           %UpdatePlot(handles)
           this.viewingWidgetHandle.plotContour=get(hObject,'Value');
       end
       % --- Executes on slider movement.
       function sliderSlice_Callback(this,hObject, ~)
           % hObject    handle to sliderSlice (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           
           % Hints: get(hObject,'Value') returns position of slider
           %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
           %UpdatePlot(handles)
           
           this.viewingWidgetHandle.slice= round(get(hObject,'Value'));
       end
       
       function radiobtnCT_Callback(this,hObject, ~)
           % hObject    handle to radiobtnCT (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           
           % Hint: get(hObject,'Value') returns toggle state of radiobtnCT
           %UpdatePlot(handles)
           this.viewingWidgetHandle.plotCT=get(hObject,'Value');
       end
       
       % --- Executes on button press in radiobtnPlan.
       function radiobtnPlan_Callback(this,hObject, ~)
           % hObject    handle to radiobtnPlan (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           
           % Hint: get(hObject,'Value') returns toggle state of radiobtnPlan
           %UpdatePlot(handles)
           this.viewingWidgetHandle.plotPlan=get(hObject,'Value');
       end
       
       % --- Executes on button press in radiobtnIsoDoseLines.
       function radiobtnIsoDoseLines_Callback(this,hObject, ~)
           % hObject    handle to radiobtnIsoDoseLines (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           
           % Hint: get(hObject,'Value') returns toggle state of radiobtnIsoDoseLines
           %UpdatePlot(handles)
           handles=this.handles;
           this.viewingWidgetHandle.plotIsoDoseLines=get(hObject,'Value');
%            if get(hObject,'Value')
%                set(handles.radiobtnIsoDoseLinesLabels,'Enable','on');
%            else
%                set(handles.radiobtnIsoDoseLinesLabels,'Enable','off');
%            end
           this.handles=handles;
       end
       
       % --- Executes on button press in radiobtnDose.
       function radiobtnDose_Callback(this,hObject, ~)
           % hObject    handle to radiobtnDose (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           
           % Hint: get(hObject,'Value') returns toggle state of radiobtnDose
           %UpdatePlot(handles)
           this.viewingWidgetHandle.plotDose=get(hObject,'Value');
       end
       
       % radio button: plot isolines labels
       function radiobtnIsoDoseLinesLabels_Callback(this,hObject, ~)
           %UpdatePlot(handles);
           this.viewingWidgetHandle.plotIsoDoseLinesLabels=get(hObject,'Value');
       end
       
       % --- Executes on button press in radioBtnIsoCenter.
       function radioBtnIsoCenter_Callback(this,hObject, eventdata)
           % hObject    handle to radioBtnIsoCenter (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           %UpdatePlot(handles)
           % Hint: get(hObject,'Value') returns toggle state of radioBtnIsoCenter
           this.viewingWidgetHandle.plotIsoCenter=get(hObject,'Value');
       end
       
%        function updateIsodoseLine(this)
%            lockState=this.viewingWidgetHandle.lockUpdate;
%            this.viewingWidgetHandle.lockUpdate=true;
%            if evalin('base','exist(''resultGUI'')')
%                 
%                 resultGUI = evalin('base','resultGUI');
%                 % select first cube if selected option does not exist
%                 if ~isfield(resultGUI,get(this.handles.popupDisplayOption,'Value'))%.viewingWidgetHandle.SelectedDisplayOption)
%                     CubeNames = fieldnames(resultGUI);
%                     dose = resultGUI.(CubeNames{1,1});
%                 else
%                     dose = resultGUI.get(this.handles.popupDisplayOption,'Value'); %(this.viewingWidgetHandle.SelectedDisplayOption);
%                 end
%                 
%                 %if function is called for the first time then set display parameters
%                 if isempty(this.viewingWidgetHandle.dispWindow{3,2})
%                     this.viewingWidgetHandle.dispWindow{3,1} = [min(dose(:)) max(dose(:))]; % set default dose range
%                     this.viewingWidgetHandle.dispWindow{3,2} = [min(dose(:)) max(dose(:))]; % set min max values
%                 end
%                 
%                 minMaxRange = this.viewingWidgetHandle.dispWindow{3,1};
%                 % if upper colorrange is defined then use it otherwise 120% iso dose
%                 upperMargin = 1;
%                 if abs((max(dose(:)) - this.viewingWidgetHandle.dispWindow{3,1}(1,2))) < 0.01  * max(dose(:))
%                     upperMargin = 1.2;
%                 end
%                 
%                 %if (length(this.viewingWidgetHandle.IsoDose_Levels) == 1 && this.viewingWidgetHandle.IsoDose_Levels(1,1) == 0) || ~this.NewIsoDoseFlag
%                     vLevels                  = [0.1:0.1:0.9 0.95:0.05:upperMargin];
%                     referenceDose            = (minMaxRange(1,2))/(upperMargin);
%                     this.viewingWidgetHandle.IsoDose_Levels   = minMaxRange(1,1) + (referenceDose-minMaxRange(1,1)) * vLevels;
%                 %end 
%             end
%            this.viewingWidgetHandle.lockUpdate=lockState;
%        end
    end
end
