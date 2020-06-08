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
            %handles=this.handles;
            this.UpdateButtonsState();
            %this.handles=handles;
        end
                
%         function viewingWidgetHandle=get.viewingWidgetHandle(this)
%             viewingWidgetHandle=this.viewingWidgetHandle;
%         end
        
        function set.viewingWidgetHandle(this,value)
            handles=this.handles;
            if isa(value,'matRad_ViewingWidget')
                this.viewingWidgetHandle=value;
                
                % get the default values from the viewer widget
                %set(handles.popupDisplayOption,'Value',this.viewingWidgetHandle.SelectedDisplayOption);
%                 if strcmp(this.viewingWidgetHandle.SelectedDisplayOption, '')
%                     set(handles.popupDisplayOption,'String','no option available');
%                 else
%                     set(handles.popupDisplayOption,'String',this.viewingWidgetHandle.SelectedDisplayOption);
%                 end

      
                if strcmp(this.viewingWidgetHandle.ProfileType,'lateral')
                    set(handles.btnProfileType,'String','longitudinal');
                else
                    set(handles.btnProfileType,'String','lateral');
                end
                
                %set(handles.btnProfileType,'String',this.viewingWidgetHandle.ProfileType);
                set(handles.popupTypeOfPlot,'Value',this.viewingWidgetHandle.typeOfPlot);
                set(handles.popupPlane,'Value',this.viewingWidgetHandle.plane);
                set(handles.radiobtnContour,'Value',this.viewingWidgetHandle.plotContour);
                set(handles.radiobtnDose,'Value',this.viewingWidgetHandle.plotDose);
                set(handles.radiobtnIsoDoseLines,'Value',this.viewingWidgetHandle.plotIsoDoseLines);
                set(handles.radiobtnIsoDoseLinesLabels,'Value',this.viewingWidgetHandle.plotIsoDoseLinesLabels);
                set(handles.radioBtnIsoCenter,'Value',this.viewingWidgetHandle.plotIsoCenter);
                set(handles.radiobtnPlan,'Value',this.viewingWidgetHandle.plotPlan);
                
                
                this.handles=handles;
                this.UpdateButtonsState();
              
                
            else
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
                'Tag','sliderSlice');
            
            h39 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','Plane Selection',...
                'Style','text',...
                'Position',[0.34258853596203 0.589041095890411 0.11969696969697 0.116438356164384],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Tag','txtPlanSelection');
            
            h40 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Slice Selection',...
                'Style','text',...
                'Position',[0.00909090909090909 0.808219178082192 0.113636363636364 0.116438356164384],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
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
                'Tag','radiobtnIsoDoseLines');
            
            h44 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Type of plot',...
                'Style','text',...
                'Position',[0.343053173241852 0.793103448275862 0.108061749571183 0.124137931034483],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
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
                'FontWeight','bold');
            
            h46 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Display option',...
                'Style','text',...
                'Position',[0.34258853596203 0.39041095890411 0.128787878787879 0.102739726027397],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
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
                'FontWeight','bold');
            
            h48 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Beam Selection',...
                'Style','text',...
                'Position',[0.00857632933104631 0.503448275862069 0.118353344768439 0.186206896551724],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
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
                'FontWeight','bold' );
            
            h51 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','GoTo',...
                'Style','text',...
                'Position',[0.604795511376553 0.801369863013698 0.0454545454545454 0.123287671232877],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Children',[],...
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
                'FontWeight','bold');
            
            h53 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot isolines labels',...
                'Style','radiobutton',...
                'Position',[0.780445969125214 0.343557168784029 0.202401372212693 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnIsoDoseLinesLabels_Callback(this,hObject,eventdata),...
                'Tag','radiobtnIsoDoseLinesLabels');
            
            h54 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'String','Offset',...
                'Style','text',...
                'Position',[0.00909090909090909 0.287104393008975 0.118181818181818 0.123287671232877],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
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
                'Tag','sliderOffset');
            
            
            h56 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','plot iso center',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.780445969125214 0.205989110707804 0.169811320754717 0.117241379310345],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radioBtnIsoCenter_Callback(this,hObject,eventdata),...
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
                'Tag','radiobtnCT');
            
            
            h59 = uicontrol(...
                'Parent',h36,...
                'Units','normalized',...
                'String','visualize plan / beams',...
                'Style','radiobutton',...
                'Value',1,...
                'Position',[0.78021978021978 0.0736842105263158 0.2002442002442 0.115789473684211],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
                'Callback',@(hObject,eventdata) radiobtnPlan_Callback(this,hObject,eventdata),...
                'Tag','radiobtnPlan');
            
            
            this.createHandles();
        end
    end
    
    methods (Access = protected)
        
        function UpdateButtonsState(this)
            handles=this.handles;
            
            % update the sliders and enable 3D view when there is data
            if evalin('base','exist(''pln'')') && evalin('base','exist(''ct'')') && ...
                    evalin('base','exist(''cst'')') && ~isempty(this.viewingWidgetHandle) %if handles.State > 0
                set(handles.btn3Dview,'Enable','on');
                 
                % set slice slider limit
                ct = evalin('base','ct');
                set(handles.sliderSlice,'Min',1,'Max',ct.cubeDim(this.viewingWidgetHandle.plane),...
                    'Value', this.viewingWidgetHandle.slice, ... %round(ct.cubeDim(this.viewingWidgetHandle.plane)/2),...
                    'SliderStep',[1/(ct.cubeDim(this.viewingWidgetHandle.plane)-1) 1/(ct.cubeDim(this.viewingWidgetHandle.plane)-1)]);
                
                pln = evalin('base','pln');
                
                if numel(pln.propStf.gantryAngles) > 1
                    % set beam selection slider
                    set(handles.sliderBeamSelection,'Min',this.viewingWidgetHandle.selectedBeam,'Max',pln.propStf.numOfBeams,...
                        'Value',this.viewingWidgetHandle.selectedBeam,...
                        'SliderStep',[1/(pln.propStf.numOfBeams-1) 1/(pln.propStf.numOfBeams-1)]);
                end
                
                % set profile offset slider
                vMinMax = [-100 100];
                vRange = sum(abs(vMinMax));
                
                ct = evalin('base','ct');
                if strcmp(get(handles.btnProfileType,'String'),'lateral')
                    SliderStep = vRange/ct.resolution.x;
                else
                    SliderStep = vRange/ct.resolution.y;
                end
                
                set(handles.sliderOffset,'Min',vMinMax(1),'Max',vMinMax(2),...
                    'Value',this.viewingWidgetHandle.profileOffset,...
                    'SliderStep',[1/SliderStep 1/SliderStep]);
                
                              
                if evalin('base','exist(''resultGUI'')')
                    set(handles.btnDVH,'Enable','on');
                    
                    Result = evalin('base','resultGUI');
                    
                    % populate the display option popup
                    if ~isempty(Result) && ~isempty(ct.cubeHU) %&& ~isfield(handles,'DispInfo')
                        
                        DispInfo = fieldnames(Result);
                        
                        for i = 1:size(DispInfo,1)
                            
                            % delete weight vectors in Result struct for plotting
                            if isstruct(Result.(DispInfo{i,1})) || isvector(Result.(DispInfo{i,1}))
                                Result = rmfield(Result,DispInfo{i,1});
                                DispInfo{i,2}=false;
                            else
                                %second dimension indicates if it should be plotted
                                DispInfo{i,2} = true;
                                % determine units
                                if strfind(DispInfo{i,1},'physicalDose')
                                    DispInfo{i,3} = '[Gy]';
                                elseif strfind(DispInfo{i,1},'alpha')
                                    DispInfo{i,3} = '[Gy^{-1}]';
                                elseif strfind(DispInfo{i,1},'beta')
                                    DispInfo{i,3} = '[Gy^{-2}]';
                                elseif strfind(DispInfo{i,1},'RBExD')
                                    DispInfo{i,3} = '[Gy(RBE)]';
                                elseif strfind(DispInfo{i,1},'LET')
                                    DispInfo{i,3} = '[keV/um]';
                                else
                                    DispInfo{i,3} = '[a.u.]';
                                end
                                DispInfo{i,4} = [];    % optional for the future: color range for plotting
                                DispInfo{i,5} = [];    % optional for the future: min max values
                            end
                        end
                        
                        set(handles.popupDisplayOption,'String',fieldnames(Result));
%                         
%                         if sum(strcmp(this.viewingWidgetHandle.SelectedDisplayOption,fieldnames(Result))) == 0
%                             set(handles.popupDisplayOption,'Value', find([DispInfo{:,2}],1,'first'));
%                             this.updateIsodoseLine();
%                             content=get(handles.popupDisplayOption,'String');
%                             this.viewingWidgetHandle.SelectedDisplayOption= content{get(handles.popupDisplayOption,'Value'),1};
%                         else
%                             set(handles.popupDisplayOption,'Value',find(strcmp(this.viewingWidgetHandle.SelectedDisplayOption,fieldnames(Result))));
%                             this.updateIsodoseLine();
%                         end

                    end
               
                else
                    set(handles.btnDVH,'Enable','off');
                    if strcmp(this.viewingWidgetHandle.SelectedDisplayOption,'')
                        set(handles.popupDisplayOption,'String','no option available');
                    else
                        set(handles.popupDisplayOption,'String',this.viewingWidgetHandle.SelectedDisplayOption);
                    end
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
                    
                    if evalin('base','exist(''pln'')')
                        pln = evalin('base','pln');
                        if numel(pln.propStf.gantryAngles) > 1  %if length(parseStringAsNum(get(handles.editGantryAngle,'String'),true)) > 1
                            
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
                    
                end
            else
                % reset slider when no data is loaded
                set(handles.sliderSlice,'Min',0,'Max',1,'Value',0,'SliderStep',[1 1]);
                % disable 3D and DVH button
                set(handles.btn3Dview,'Enable','off');
                set(handles.btnDVH,'Enable','off');
                set(handles.popupDisplayOption,'String','no option available');
                
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
                
                if evalin('base','exist(''pln'')') && evalin('base','exist(''ct'')') && ...
                        evalin('base','exist(''cst'')') %if handles.State > 0
                    pln = evalin('base','pln');
                    if numel(pln.propStf.gantryAngles) > 1  %if length(parseStringAsNum(get(handles.editGantryAngle,'String'),true)) > 1
                        
                        set(handles.sliderBeamSelection,'Enable','on');
                        %this.viewingWidgetHandle.selectedBeam = 1;
                        %pln = evalin('base','pln');
%                         set(handles.sliderBeamSelection,'Min',this.viewingWidgetHandle.selectedBeam,'Max',pln.propStf.numOfBeams,...
%                             'Value',this.viewingWidgetHandle.selectedBeam,...
%                             'SliderStep',[1/(pln.propStf.numOfBeams-1) 1/(pln.propStf.numOfBeams-1)],...
%                             'Enable','on');
%                         
%                     else
%                         this.viewingWidgetHandle.selectedBeam = 1;
                    end
                    
                    %this.viewingWidgetHandle.profileOffset = get(handles.sliderOffset,'Value');
                    
%                     vMinMax = [-100 100];
%                     vRange = sum(abs(vMinMax));
%                     
%                     ct = evalin('base','ct');
%                     if strcmp(get(handles.btnProfileType,'String'),'lateral')
%                         SliderStep = vRange/ct.resolution.x;
%                     else
%                         SliderStep = vRange/ct.resolution.y;
%                     end
                    
%                     set(handles.sliderOffset,'Min',vMinMax(1),'Max',vMinMax(2),...
%                         'Value',this.viewingWidgetHandle.profileOffset,...
%                         'SliderStep',[1/SliderStep 1/SliderStep],...
%                         'Enable','on');
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
            this.updateIsodoseLine();
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
            handles = this.handles;
            
           
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
                this.handles = handles;
            end
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
            handles = this.handles;
            this.viewingWidgetHandle.profileOffset = get(hObject,'Value');
            %UpdatePlot(handles);
            this.handles = handles;
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
       
       function updateIsodoseLine(this)
           this.viewingWidgetHandle.lockUpdate=true;
           if evalin('base','exist(''resultGUI'')')
                
                resultGUI = evalin('base','resultGUI');
                % select first cube if selected option does not exist
                if ~isfield(resultGUI,get(this.handles.popupDisplayOption,'Value'))%.viewingWidgetHandle.SelectedDisplayOption)
                    CubeNames = fieldnames(resultGUI);
                    dose = resultGUI.(CubeNames{1,1});
                else
                    dose = resultGUI.get(this.handles.popupDisplayOption,'Value'); %(this.viewingWidgetHandle.SelectedDisplayOption);
                end
                
                %if function is called for the first time then set display parameters
                if isempty(this.viewingWidgetHandle.dispWindow{3,2})
                    this.viewingWidgetHandle.dispWindow{3,1} = [min(dose(:)) max(dose(:))]; % set default dose range
                    this.viewingWidgetHandle.dispWindow{3,2} = [min(dose(:)) max(dose(:))]; % set min max values
                end
                
                minMaxRange = this.viewingWidgetHandle.dispWindow{3,1};
                % if upper colorrange is defined then use it otherwise 120% iso dose
                upperMargin = 1;
                if abs((max(dose(:)) - this.viewingWidgetHandle.dispWindow{3,1}(1,2))) < 0.01  * max(dose(:))
                    upperMargin = 1.2;
                end
                
                %if (length(this.viewingWidgetHandle.IsoDose_Levels) == 1 && this.viewingWidgetHandle.IsoDose_Levels(1,1) == 0) || ~this.NewIsoDoseFlag
                    vLevels                  = [0.1:0.1:0.9 0.95:0.05:upperMargin];
                    referenceDose            = (minMaxRange(1,2))/(upperMargin);
                    this.viewingWidgetHandle.IsoDose_Levels   = minMaxRange(1,1) + (referenceDose-minMaxRange(1,1)) * vLevels;
                %end 
            end
           this.viewingWidgetHandle.lockUpdate=false;
       end
    end
end
