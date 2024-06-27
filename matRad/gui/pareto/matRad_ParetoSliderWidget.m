classdef matRad_ParetoSliderWidget < matRad_Widget
    % matRad_ParetoSliderWidget class to generate GUI widget to move sliders for
    % Pareto surface navigation
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
        paretoSliceWidgetHandle
        objectiveSliders
        fixationButtons
        helperObject
        optiProb
    end
    
    methods
        function this = matRad_ParetoSliderWidget(handleParent,paretoSliceWidgetHandle)

            matRad_cfg = MatRad_Config.instance();
            this = this@matRad_Widget(handleParent);

            this.paretoSliceWidgetHandle = paretoSliceWidgetHandle;
            this.objectiveSliders = {};
            %TODO: Best way to pass objective function values?
            %fAll = evalin('base','retStruct.finds');
            %wAll = evalin('base','retStruct.weights');
            %this.optiProb = evalin('base','retStruct.optiProb');
            %this.helperObject = matRad_ParetoData(fAll,wAll);
            this.initialize();
            this.plotSlice(evalin('base','ParetoHelperObject.currentWeights'));
        end
                              
    end
    
    methods (Access = protected)
        function this = createLayout(this)
            h1 = this.widgetHandle;

            matRad_cfg = MatRad_Config.instance();
            % handle environment
            switch matRad_cfg.env
                case 'MATLAB'
                    set(h1,'SizeChangedFcn',@(hObject,eventdata) widget_SizeChangedFcn(this,hObject,eventdata));
                case 'OCTAVE'
                    % this creates an infinite loop in octave
                    %set(h1,'SizeChangedFcn',@(hObject,eventdata) widget_SizeChangedFcn(this,hObject,eventdata));
            end

            this.createHandles();
            
        end

        function this = doUpdate(this,evt)
            doUpdate = true;

            if nargin == 2
                %At pln changes and at cst/cst (for Isocenter and new settings) 
                %we need to update
                doUpdate = this.checkUpdateNecessary({'cst'},evt);
            end
            
            if doUpdate
                if evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')
                    generateSliderTable(this, evalin('base','cst'));
                else
                    delete(get(this.widgetHandle,'Children'));
                end
            end
           
        end
    end
    
    methods(Access = private)
        
        function cst = generateSliderTable(this,cst)
            matRad_cfg = MatRad_Config.instance();
            ParetoHelperObject = evalin('base','ParetoHelperObject');
            %cst = updateStructureTable(this,cst);
            handles = this.handles;
            SliderPanel = this.widgetHandle;
            
            SliderPanelPos = get(SliderPanel,'Position');
            SliderPanelPosUnit = get(SliderPanel,'Units');
            set(SliderPanel,'Units','pixels');
            SliderPanelPosPix = get(SliderPanel,'Position');  
            set(SliderPanel,'Units',SliderPanelPosUnit);
            aspectRatio = SliderPanelPosPix(3) / SliderPanelPosPix(4);
            
            %Parameters for line height
            objHeight = 0.055;% 22;
            lineHeight = 0.1; %25; %Height of a table line
            yTopSep = 0.12;%40; %Separation of the first line from the top
            %tableViewHeight = SliderPanelPos(4) - yTopSep; %Full height of the view
            tableViewHeight = 1 - yTopSep;
            
            %Widths of the fields
            buttonW = objHeight*1.5 / aspectRatio; % Make button squared
            nameW = 3*buttonW;%60;
            sliderW = 5*buttonW;
            functionW = 5*buttonW;%120;
            paramTitleW = 4*buttonW;%120;
            paramW = 1*buttonW;%30;
            fieldSep = 0.25*buttonW; %Separation between fields horizontally

            %disp(num2str(sliderPos));
            SliderPanelChildren = get(SliderPanel,'Children');
            cstVertTableScroll = findobj(SliderPanelChildren,'Tag','VerticalSlider');
            if isempty(cstVertTableScroll)
                sliderPos = 0;
            else
                sliderPos = get(cstVertTableScroll,'Max') - get(cstVertTableScroll,'Value');
            end
            
            ypos = @(c) tableViewHeight - c*lineHeight + sliderPos;
            
            delete(SliderPanelChildren);


            %Creates a dummy axis to allow for the use of textboxes instead of uicontrol to be able to use the (la)tex interpreter
            tmpAxes = axes('Parent',SliderPanel,'units','normalized','position',[0 0 1 1],'visible','off', 'FontSize',8);
            
            organTypes = {'OAR', 'TARGET','IGNORED'};
            
            %Get all Classes & classNames
            classNames = matRad_getObjectivesAndConstraints();
            
            numOfObjectives = 0;
            for i = 1:size(cst,1)
                for j = 1:numel(cst{i,6})
                    if isa(cst{i,6}{j},'DoseObjectives.matRad_DoseObjective')
                        numOfObjectives = numOfObjectives + 1;
                    end
                end
            end
            %line 
            % step  count 
            cnt = 0;
            
            newline = '\n';
            
            %Setup Headlines
            xPos = 0.01; %5

            h2 = uicontrol(SliderPanel,'Style','text', ...
                'String','VOI name', ...
                'Units','normalized', ...
                'Position',[xPos ypos(cnt) nameW objHeight], ...
                'FontSize',matRad_cfg.gui.fontSize, ...
                'Tooltip','Name of the structure with objective/constraint', ...
                'BackgroundColor',matRad_cfg.gui.backgroundColor, ...
                'ForegroundColor',matRad_cfg.gui.textColor);

            tmp_pos = get(h2,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;

            h7 = uicontrol(SliderPanel,'Style','text', ...
                'String','Function slider', ...
                'Units','normalized', ...
                'Position',[xPos ypos(cnt) functionW objHeight], ...
                'FontSize',matRad_cfg.gui.fontSize, ...
                'Tooltip','Objective/Constraint function type', ...
                'BackgroundColor',matRad_cfg.gui.backgroundColor, ...
                'ForegroundColor',matRad_cfg.gui.textColor);

            tmp_pos = get(h7,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;

            h3 = uicontrol(SliderPanel,'Style','text', ...
                'String','Fix slider', ...
                'Units','normalized', ...
                'Position',[xPos ypos(cnt) buttonW objHeight], ...
                'FontSize',matRad_cfg.gui.fontSize, ...
                'Tooltip','Button to fixate a slider', ...
                'BackgroundColor',matRad_cfg.gui.backgroundColor, ...
                'ForegroundColor',matRad_cfg.gui.textColor);

            tmp_pos = get(h3,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            
            h6 = uicontrol(SliderPanel,'Style','text', ...
                'String','Function', ...
                'Units','normalized', ...
                'Position',[xPos ypos(cnt) sliderW objHeight], ...
                'FontSize',matRad_cfg.gui.fontSize, ...
                'Tooltip','Objective/Constraint function type', ...
                'BackgroundColor',matRad_cfg.gui.backgroundColor, ...
                'ForegroundColor',matRad_cfg.gui.textColor);

            tmp_pos = get(h6,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;

            h8 = uicontrol(SliderPanel,'Style','text', ...
                'String','| Parameters', ...
                'Units','normalized', ...
                'Position',[xPos ypos(cnt) paramTitleW objHeight], ...
                'FontSize',matRad_cfg.gui.fontSize, ...
                'Tooltip','List of parameters', ...
                'HorizontalAlignment','left', ...
                'BackgroundColor',matRad_cfg.gui.backgroundColor, ...
                'ForegroundColor',matRad_cfg.gui.textColor);

            tmp_pos = get(h8,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            cnt = cnt + 1;
            %}
            %Create Objectives / Constraints controls
            for i = 1:size(cst,1)
                if strcmp(cst(i,3),'IGNORED')~=1
                    %Compatibility Layer for old objective format
                    if isstruct(cst{i,6})
                        cst{i,6} = num2cell(arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6}));
                    end
                    for j=1:numel(cst{i,6})
                        %TODO: Add check if constraint or objective
                        obj = cst{i,6}{j};
                        
                        %Convert to class if not
                        if ~isa(obj,'matRad_DoseOptimizationFunction')
                            try
                                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                            catch ME
                                this.showWarning('Objective/Constraint not valid!\n%s',ME.message)
                                continue;
                            end
                        end
                        
                        if isa(obj,'DoseObjectives.matRad_DoseObjective')
                            xPos = 0.01;
                            h = uicontrol(SliderPanel','Style','edit', ...
                                'String',cst{i,2}, ...
                                'Units','normalized', ...
                                'Position',[xPos ypos(cnt) nameW objHeight], ...
                                'FontSize',matRad_cfg.gui.fontSize, ...
                                'FontName',matRad_cfg.gui.fontName, ...
                                'FontWeight',matRad_cfg.gui.fontWeight, ...
                                'BackgroundColor',matRad_cfg.gui.elementColor, ...
                                'ForegroundColor',matRad_cfg.gui.textColor, ...
                                'Tooltip','Name',...
                                'Enable','inactive',... %Disable editing of name atm
                                'UserData',[i,2]);
                            
                            tmp_pos = get(h,'Position');
                            xPos = xPos + tmp_pos(3) + fieldSep;
                            %Slider
                            this.objectiveSliders{cnt} = uicontrol(SliderPanel,'Units','normalized',...
                                'String','Slider',...
                                'Tooltip','Adjust objective penalty',...
                                'Position',[xPos ypos(cnt) sliderW objHeight], ...
                                'Style','slider',...
                                'BusyAction','cancel',...
                                'Interruptible','off',...
                                'BackgroundColor',matRad_cfg.gui.elementColor,...
                                'ForegroundColor',matRad_cfg.gui.textColor,...
                                'FontSize',matRad_cfg.gui.fontSize,...
                                'FontName',matRad_cfg.gui.fontName,...
                                'FontWeight',matRad_cfg.gui.fontWeight, ...
                                'Callback',@(hObject,eventdata)ObjectiveSlider_callback(this,hObject,eventdata,cnt));
                            
                            this.objectiveSliders{cnt}.Value = ParetoHelperObject.currentPoint(cnt);
                
                            tmp_pos = get(this.objectiveSliders{cnt},'Position');
                            xPos = xPos + tmp_pos(3) + fieldSep;

                            this.fixationButtons{cnt} = uicontrol(SliderPanel,'Style','pushbutton', ...
                            'String','fix', ...
                            'Units','normalized', ...
                            'Position',[xPos ypos(cnt) buttonW objHeight], ...
                            'FontSize',matRad_cfg.gui.fontSize, ...
                            'FontName',matRad_cfg.gui.fontName, ...
                            'FontWeight',matRad_cfg.gui.fontWeight, ...
                            'BackgroundColor',matRad_cfg.gui.elementColor, ...
                            'ForegroundColor',matRad_cfg.gui.textColor, ...
                            'Tooltip','Remove Objective/Constraint',...
                            'Callback',@(hObject,eventdata)FixButton_callback(this,hObject,eventdata,cnt));
                            tmp_pos = get(this.fixationButtons{cnt},'Position');
                            xPos = xPos + tmp_pos(3) + fieldSep;

                            h = uicontrol(SliderPanel,'Style','text', ...
                                'String',obj.name, ...
                                'Units','normalized', ...
                                'Position',[xPos ypos(cnt) functionW objHeight], ...
                                'FontSize',matRad_cfg.gui.fontSize, ...
                                'FontName',matRad_cfg.gui.fontName, ...
                                'FontWeight',matRad_cfg.gui.fontWeight, ...
                                'BackgroundColor',matRad_cfg.gui.elementColor, ...
                                'ForegroundColor',matRad_cfg.gui.textColor, ...
                                'Tooltip','Objective/Constraint',...
                                'UserData',{[i,j],classNames(1,:)});
    
                            tmp_pos = get(h,'Position');
                            xPos = xPos + tmp_pos(3) + fieldSep;
                            
                           
                            
                            for p = 1:numel(obj.parameterNames)
    %                             h = text('Parent',tmpAxes,'String',['| ' obj.parameterNames{p} ':'],'VerticalAlignment','middle','Units','normalized','Position',[xPos ypos(cnt)+lineHeight/2],'Interpreter','tex','FontWeight','normal',...
    %                                 'FontSize',get(SliderPanel,'FontSize'),'FontName',get(SliderPanel,'FontName'),'FontUnits',get(SliderPanel,'FontUnits'),'FontWeight','normal');%[xPos ypos(cnt) 100 objHeight]);
                                % there is no fontsize for SliderPanel
                                h = text('Parent',tmpAxes, ...
                                    'String',['| ' obj.parameterNames{p} ':'], ...
                                    'VerticalAlignment','middle', ...
                                    'Units','normalized', ...
                                    'Position',[xPos ypos(cnt)+lineHeight/2], ...
                                    'Interpreter','tex', ...
                                    'FontSize',matRad_cfg.gui.fontSize, ...
                                    'FontName',matRad_cfg.gui.fontName, ...
                                    'FontWeight',matRad_cfg.gui.fontWeight, ...
                                    'BackgroundColor',matRad_cfg.gui.backgroundColor, ...
                                    'Color',matRad_cfg.gui.textColor);%[xPos ypos(cnt) 100 objHeight]);
    
                                tmp_pos = get(h,'Extent');
                                xPos = xPos + tmp_pos(3) + fieldSep;
                                %h = annotation(SliderPanel,'textbox','String',obj.parameters{1,p},'Units','pix','Position', [xPos ypos(cnt) 100 objHeight],'Interpreter','Tex');
                                
                                %Check if we have a cell and therefore a parameter list
                                if iscell(obj.parameterTypes{p})
                                    h = uicontrol(SliderPanel,'Style','popupmenu', ...
                                        'String',obj.parameterTypes{p}', ...
                                        'Value',obj.parameters{p}, ...
                                        'Tooltip',obj.parameterNames{p}, ...
                                        'Units','normalized', ...
                                        'Position',[xPos ypos(cnt) paramW*2 objHeight], ...
                                        'FontSize',matRad_cfg.gui.fontSize, ...
                                        'FontName',matRad_cfg.gui.fontName, ...
                                        'FontWeight',matRad_cfg.gui.fontWeight, ...
                                        'BackgroundColor',matRad_cfg.gui.elementColor, ...
                                        'ForegroundColor',matRad_cfg.gui.textColor, ...
                                        'UserData',[i,j,p]);
                                else
                                    h = uicontrol(SliderPanel,'Style','text', ...
                                        'String',num2str(obj.parameters{p}), ...
                                        'Tag','FunctionValueSlider',...
                                        'Tooltip',obj.parameterNames{p}, ...
                                        'Units','normalized', ...
                                        'Position',[xPos ypos(cnt) paramW objHeight], ...
                                        'FontSize',matRad_cfg.gui.fontSize, ...
                                        'FontName',matRad_cfg.gui.fontName, ...
                                        'FontWeight',matRad_cfg.gui.fontWeight, ...
                                        'BackgroundColor',matRad_cfg.gui.elementColor, ...
                                        'ForegroundColor',matRad_cfg.gui.textColor, ...
                                        'UserData',[i,j,p]);
                                end
                                
                                tmp_pos = get(h,'Position');
                                xPos = xPos + tmp_pos(3) + fieldSep;
                            end


                            cnt = cnt +1;
                        end     
                    end
                end
            end
            %Calculate Scrollbar
            lastPos = ypos(cnt);
            firstPos = ypos(0);
            tableHeight = abs(firstPos - lastPos);
            
            exceedFac = tableHeight / tableViewHeight;
            if exceedFac > 1
                sliderFac = exceedFac - 1;
                uicontrol(SliderPanel,'Style','slider', ...
                    'Units','normalized', ...
                    'Tag','VerticalSlider',...
                    'Position',[0.975 0 0.025 1], ...
                    'FontSize',matRad_cfg.gui.fontSize, ...
                    'FontName',matRad_cfg.gui.fontName, ...
                    'FontWeight',matRad_cfg.gui.fontWeight, ...
                    'BackgroundColor',matRad_cfg.gui.elementColor, ...
                    'ForegroundColor',matRad_cfg.gui.textColor, ...
                    'Min',0,'Max',ceil(sliderFac)*tableViewHeight, ...
                    'SliderStep',[lineHeight tableViewHeight] ./ (ceil(sliderFac)*tableViewHeight), ...
                    'Value',ceil(sliderFac)*tableViewHeight - sliderPos, ...
                    'Callback', @(hObject,eventdata)cstTableSlider_Callback(this,hObject,eventdata));
            end
            this.handles = handles;
        end        

        function ObjectiveSlider_callback(this,slider,~,idx)
            % 
            ParetoHelperObject = evalin('base','ParetoHelperObject');
            v = matRad_paretoSurfaceNavigation(ParetoHelperObject.allPoints',ParetoHelperObject.currentPoint,slider.Value,idx,ParetoHelperObject.availableUpbound);
            
            if numel(v) == 0
            %navigation algorithm didnt find a new point
                for i = 1:numel(this.objectiveSliders)
                    set(this.objectiveSliders{i},'Value',ParetoHelperObject.currentPoint(i));
                end
            else
                wNew = ParetoHelperObject.allWeights*v;

                fNew = ParetoHelperObject.allPoints'*v;% evaluate function(strictly necessary?)
                %fNew = this.optiProb.normalizeObjectives(fNew');
                

                ParetoHelperObject.currentPoint = fNew;
                ParetoHelperObject.currentWeights = wNew;
                
                %Open: Current
                for i = 1:numel(this.objectiveSliders)
                    set(this.objectiveSliders{i},'Value',fNew(i));
                end
                
                this.plotSlice(wNew);
            end
        end

        function plotSlice(this,w)
            slice = round(evalin('base','pln.propStf.isoCenter(1,3)./ct.resolution.z'));
            cubes = matRad_calcFastCubes(w,evalin('base','dij'),evalin('base','pln'));
            
            matRad_plotSliceWrapper(this.paretoSliceWidgetHandle.DosePlotAxes,evalin('base','ct'),evalin('base','cst'),1,cubes,3,slice,[],[],[],[],[],[],[],[],[],'LineWidth',2);
        end

        function FixButton_callback(this,~,~,idx)
            
            ParetoHelperObject = evalin('base','ParetoHelperObject');
            [lb,ub] = ParetoHelperObject.restrictObjective(idx,this.objectiveSliders{idx}.Value); %update refObjects bounds
            
    
            for i = 1:numel(this.objectiveSliders)
                set(this.objectiveSliders{i},'Min',lb(i));
                set(this.objectiveSliders{i},'Max',ub(i));
            end
        end
    

                % --- Executes when the widget is resized.
        function widget_SizeChangedFcn(this,hObject, eventdata)
            % hObject    handle to h1 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)      
            try
                generateSliderTable(this,evalin('base','cst'));
            catch
            end
        end

                % --- Executes on slider movement.
        function cstTableSlider_Callback(this,~,~)
            % hObject    handle to cstTableSlider (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: get(hObject,'Value') returns position of slider
            %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
            try
                generateSliderTable(this,evalin('base','cst'));
            catch
            end
        end
    end
end
