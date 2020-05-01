classdef matRad_OptimizationWidget < matRad_Widget
    
    properties
        
    end
    
    methods
        function this = matRad_OptimizationWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[100 45 90 15],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Optimization',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1');
                    %'PaperSize',[20.99999864 29.69999902]
                
                
            end
            this = this@matRad_Widget(handleParent);
        end
       
        function this=initialize(this)
             this.update();
            
        end
        
        function this=update(this)
            
            if evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')
                generateCstTable(this, evalin('base','cst'));
            else
            end
           
        end
        
        
        
        function changeWorkspace(obj)
            notify(obj, 'workspaceChanged');
        end
    end
    
    methods (Access = protected)
        function this = createLayout(this)
            h1 = this.widgetHandle;
            set(h1,'SizeChangedFcn',@(hObject,eventdata) widget_SizeChangedFcn(this,hObject,eventdata));
      
            this.createHandles();
            
        end
    end
    
    methods(Access = private)
        
        function cst = generateCstTable(this,cst)
            %cst = updateStructureTable(this,cst);
            handles = this.handles;
            cstPanel = this.widgetHandle;
            
            cstPanelPos = get(cstPanel,'Position');
            cstPanelPosUnit = get(cstPanel,'Units');
            set(cstPanel,'Units','pixels');
            cstPanelPosPix = get(cstPanel,'Position');
            set(cstPanel,'Units',cstPanelPosUnit);
            aspectRatio = cstPanelPosPix(3) / cstPanelPosPix(4);
            
            %Parameters for line height
            objHeight = 0.095;% 22;
            lineHeight = 0.1; %25; %Height of a table line
            yTopSep = 0.12;%40; %Separation of the first line from the top
            %tableViewHeight = cstPanelPos(4) - yTopSep; %Full height of the view
            tableViewHeight = 1 - yTopSep;
            
            %Widths of the fields
            buttonW = objHeight / aspectRatio; % Make button squared
            nameW = 3.5*buttonW;%60;
            typeW = 3*buttonW;%70;
            opW = buttonW;%objHeight;
            functionW = 6*buttonW;%120;
            penaltyW = 2*buttonW;%40;
            paramTitleW = 4*buttonW;%120;
            paramW = 2*buttonW;%30;
            fieldSep = 0.25*buttonW; %Separation between fields horizontally
            
            %Scrollbar
            cstPanelChildren = get(cstPanel,'Children');
            cstVertTableScroll = findobj(cstPanelChildren,'Style','slider');
            if isempty(cstVertTableScroll)
                sliderPos = 0;
            else
                sliderPos = get(cstVertTableScroll,'Max') - get(cstVertTableScroll,'Value');
            end
            %disp(num2str(sliderPos));
            ypos = @(c) tableViewHeight - c*lineHeight + sliderPos;
            
            delete(cstPanelChildren);
            
            %Creates a dummy axis to allow for the use of textboxes instead of uicontrol to be able to use the (la)tex interpreter
            tmpAxes = axes('Parent',cstPanel,'units','normalized','position',[0 0 1 1],'visible','off');
            
            organTypes = {'OAR', 'TARGET'};
            
            %Get all Classes & classNames
            classNames = matRad_getObjectivesAndConstraints();
            
            numOfObjectives = 0;
            for i = 1:size(cst,1)
                if ~isempty(cst{i,6})
                    numOfObjectives = numOfObjectives + numel(cst{i,6});
                end
            end
            
            cnt = 0;
            
            newline = '\n';
            
            %Setup Headlines
            xPos = 0.01; %5
            
            h = uicontrol(cstPanel,'Style','text','String','+/-','Units','normalized','Position',[xPos ypos(cnt) buttonW objHeight],'TooltipString','Remove or add Constraint or Objective');
            tmp_pos = get(h,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            h = uicontrol(cstPanel,'Style','text','String','VOI name','Units','normalized','Position',[xPos ypos(cnt) nameW objHeight],'TooltipString','Name of the structure with objective/constraint');
            tmp_pos = get(h,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            h = uicontrol(cstPanel,'Style','text','String','VOI type','Units','normalized','Position',[xPos ypos(cnt) typeW objHeight],'TooltipString','Segmentation Classification');
            tmp_pos = get(h,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            h = uicontrol(cstPanel,'Style','text','String','OP','Units','normalized','Position',[xPos ypos(cnt) opW objHeight],'TooltipString',['Overlap Priority' char(10) '(Smaller number overlaps higher number)']);
            tmp_pos = get(h,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            h = uicontrol(cstPanel,'Style','text','String','Function','Units','normalized','Position',[xPos ypos(cnt) functionW objHeight],'TooltipString','Objective/Constraint function type');
            tmp_pos = get(h,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            h = uicontrol(cstPanel,'Style','text','String','p','Units','normalized','Position',[xPos ypos(cnt) penaltyW objHeight],'TooltipString','Optimization penalty');
            tmp_pos = get(h,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            h = uicontrol(cstPanel,'Style','text','String','| Parameters','Units','normalized','Position',[xPos ypos(cnt) paramTitleW objHeight],'TooltipString','List of parameters','HorizontalAlignment','left');
            tmp_pos = get(h,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            cnt = cnt + 1;
            
            %Create Objectives / Constraints controls
            for i = 1:size(cst,1)
                if strcmp(cst(i,3),'IGNORED')~=1
                    %Compatibility Layer for old objective format
                    if isstruct(cst{i,6})
                        cst{i,6} = num2cell(arrayfun(@matRad_DoseOptimizationFunction.convertOldOptimizationStruct,cst{i,6}));
                    end
                    for j=1:numel(cst{i,6})
                        
                        obj = cst{i,6}{j};
                        
                        %Convert to class if not
                        if ~isa(obj,'matRad_DoseOptimizationFunction')
                            try
                                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                            catch ME
                                warning('Objective/Constraint not valid!\n%s',ME.message)
                                continue;
                            end
                        end
                        
                        xPos = 0.01;%5;
                        
                        h = uicontrol(cstPanel,'Style','pushbutton','String','-','Units','normalized','Position',[xPos ypos(cnt) buttonW objHeight],'TooltipString','Remove Objective/Constraint','Callback',@(hObject,eventdata)btObjRemove_Callback(this,hObject,eventdata),...
                            'UserData',[i,j]);
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        h = uicontrol(cstPanel','Style','edit','String',cst{i,2},'Units','normalized','Position',[xPos ypos(cnt) nameW objHeight],'TooltipString','Name',...
                            'Enable','inactive',... %Disable editing of name atm
                            'UserData',[i,2],'Callback',@(hObject,eventdata)editCstParams_Callback(this,hObject,eventdata)); %Callback added, however, editing is disabled atm
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        h = uicontrol(cstPanel,'Style','popupmenu','String',organTypes','Value',find(strcmp(cst{i,3},organTypes)),'Units','normalized','Position',[xPos ypos(cnt) typeW objHeight],'TooltipString','Segmentation Classification',...
                            'UserData',[i,3],'Callback',@(hObject,eventdata)editCstParams_Callback(this,hObject,eventdata));
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        h = uicontrol(cstPanel,'Style','edit','String',num2str(cst{i,5}.Priority),'Units','normalized','Position',[xPos ypos(cnt) opW objHeight],'TooltipString',['Overlap Priority' newline '(Smaller number overlaps higher number)'],...
                            'UserData',[i,5],'Callback',@(hObject,eventdata)editCstParams_Callback(this,hObject,eventdata));
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        
                        h = uicontrol(cstPanel,'Style','popupmenu','String',classNames(2,:)','Value',find(strcmp(obj.name,classNames(2,:))),'Units','normalized','Position',[xPos ypos(cnt) functionW objHeight],'TooltipString','Select Objective/Constraint',...
                            'UserData',{[i,j],classNames(1,:)},'Callback',@(hObject,eventdata)changeObjFunction_Callback(this,hObject,eventdata));
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        
                        %Check if we have an objective to display penalty
                        if isa(obj,'DoseObjectives.matRad_DoseObjective')
                            h = uicontrol(cstPanel,'Style','edit','String',num2str(obj.penalty),'Units','normalized','Position',[xPos ypos(cnt) penaltyW objHeight],'TooltipString','Objective Penalty','UserData',[i,j,0],...
                                'Callback',@(hObject,eventdata)editObjParam_Callback(this,hObject,eventdata));
                        else
                            h = uicontrol(cstPanel,'Style','edit','String','----','Units','normalized','Position',[xPos ypos(cnt) penaltyW objHeight],'Enable','off');
                        end
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        
                        for p = 1:numel(obj.parameterNames)
%                             h = text('Parent',tmpAxes,'String',['| ' obj.parameterNames{p} ':'],'VerticalAlignment','middle','Units','normalized','Position',[xPos ypos(cnt)+lineHeight/2],'Interpreter','tex','FontWeight','normal',...
%                                 'FontSize',get(cstPanel,'FontSize'),'FontName',get(cstPanel,'FontName'),'FontUnits',get(cstPanel,'FontUnits'),'FontWeight','normal');%[xPos ypos(cnt) 100 objHeight]);
                            % there is no fontsize for cstPanel
                            h = text('Parent',tmpAxes,'String',['| ' obj.parameterNames{p} ':'],'VerticalAlignment','middle','Units','normalized','Position',[xPos ypos(cnt)+lineHeight/2],'Interpreter','tex','FontWeight','normal',...
                                'FontWeight','normal');%[xPos ypos(cnt) 100 objHeight]);
                            tmp_pos = get(h,'Extent');
                            xPos = xPos + tmp_pos(3) + fieldSep;
                            %h = annotation(cstPanel,'textbox','String',obj.parameters{1,p},'Units','pix','Position', [xPos ypos(cnt) 100 objHeight],'Interpreter','Tex');
                            
                            %Check if we have a cell and therefore a parameter list
                            if iscell(obj.parameterTypes{p})
                                h = uicontrol(cstPanel,'Style','popupmenu','String',obj.parameterTypes{p}','Value',obj.parameters{p},'TooltipString',obj.parameterNames{p},'Units','normalized','Position',[xPos ypos(cnt) paramW*2 objHeight],'UserData',[i,j,p],...
                                    'Callback',@(hObject,eventdata)editObjParam_Callback(this,hObject,eventdata));
                            else
                                h = uicontrol(cstPanel,'Style','edit','String',num2str(obj.parameters{p}),'TooltipString',obj.parameterNames{p},'Units','normalized','Position',[xPos ypos(cnt) paramW objHeight],'UserData',[i,j,p],...
                                    'Callback',@(hObject,eventdata)editObjParam_Callback(this,hObject,eventdata));
                            end
                            
                            tmp_pos = get(h,'Position');
                            xPos = xPos + tmp_pos(3) + fieldSep;
                        end
                        
                        cnt = cnt +1;
                    end
                end
            end
            xPos = 0.01; %5
            hAdd = uicontrol(cstPanel,'Style','pushbutton','String','+','Units','normalized','Position',[xPos ypos(cnt) buttonW objHeight],...
                'TooltipString','Add Objective/Constraint','Callback',@(hObject,eventdata)btObjAdd_Callback(this,hObject,eventdata)); %{@btObjAdd_Callback,handles});
            tmp_pos = get(hAdd,'Position');
            xPos = xPos + tmp_pos(3) + fieldSep;
            h = uicontrol(cstPanel,'Style','popupmenu','String',cst(:,2)','Units','normalized','Position',[xPos ypos(cnt) nameW objHeight]);
            set(hAdd,'UserData',h);
            
            %Calculate Scrollbar
            lastPos = ypos(cnt);
            firstPos = ypos(0);
            tableHeight = abs(firstPos - lastPos);
            
            exceedFac = tableHeight / tableViewHeight;
            if exceedFac > 1
                sliderFac = exceedFac - 1;
                uicontrol(cstPanel,'Style','slider','Units','normalized','Position',[0.975 0 0.025 1],'Min',0,'Max',ceil(sliderFac)*tableViewHeight,'SliderStep',[lineHeight tableViewHeight] ./ (ceil(sliderFac)*tableViewHeight),'Value',ceil(sliderFac)*tableViewHeight - sliderPos,'Callback',{@cstTableSlider_Callback,handles});
            end
            
            this.handles = handles;
        end
        
%         function cst = updateStructureTable(this,cst)
%                 
%             handles = this.handles;
%             colorAssigned = true;
%             
%             % check whether all structures have an assigned color
%             for i = 1:size(cst,1)
%                 if ~isfield(cst{i,5},'visibleColor')
%                     colorAssigned = false;
%                     break;
%                 elseif isempty(cst{i,5}.visibleColor)
%                     colorAssigned = false;
%                     break;
%                 end
%             end
%             
%             
%             % assign color if color assignment is not already present or inconsistent
%             if colorAssigned == false
%                 m = 64;
%                 colorStep = ceil(m/size(cst,1));
%                 colors    = colorcube(colorStep*size(cst,1));
%                 % spread individual VOI colors in the colorcube color palette
%                 colors    = colors(1:colorStep:end,:);
%                 
%                 for i = 1:size(cst,1)
%                     cst{i,5}.visibleColor = colors(i,:);
%                 end
%             end
%             
%             this.handles= handles;
%             changeWorkspace(this);
%         end
        
%         % precompute contours of VOIs
%         function cst = precomputeContours(this,ct,cst)
%             mask = zeros(ct.cubeDim); % create zero cube with same dimeonsions like dose cube
%             for s = 1:size(cst,1)
%                 cst{s,7} = cell(max(ct.cubeDim(:)),3);
%                 mask(:) = 0;
%                 mask(cst{s,4}{1}) = 1;
%                 for slice = 1:ct.cubeDim(1)
%                     if sum(sum(mask(slice,:,:))) > 0
%                         cst{s,7}{slice,1} = contourc(squeeze(mask(slice,:,:)),.5*[1 1]);
%                     end
%                 end
%                 for slice = 1:ct.cubeDim(2)
%                     if sum(sum(mask(:,slice,:))) > 0
%                         cst{s,7}{slice,2} = contourc(squeeze(mask(:,slice,:)),.5*[1 1]);
%                     end
%                 end
%                 for slice = 1:ct.cubeDim(3)
%                     if sum(sum(mask(:,:,slice))) > 0
%                         cst{s,7}{slice,3} = contourc(squeeze(mask(:,:,slice)),.5*[1 1]);
%                     end
%                 end
%             end
%         end

        % --- Executes when the widget is resized.
        function widget_SizeChangedFcn(this,hObject, eventdata)
            % hObject    handle to h1 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            try
                generateCstTable(this,evalin('base','cst'));
            catch
            end
        end
        
        function btObjAdd_Callback(this,hObject, ~)
            % hObject    handle to btnuiTableAdd (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles=this.handles;
            
            popupHandle = get(hObject,'UserData');
            cstIndex = get(popupHandle,'Value');
            
            cst = evalin('base','cst');
            %Add Standard Objective
            if strcmp(cst{cstIndex,3},'TARGET')
                cst{cstIndex,6}{end+1} = struct(DoseObjectives.matRad_SquaredDeviation);
            else
                cst{cstIndex,6}{end+1} = struct(DoseObjectives.matRad_SquaredOverdosing);
            end
            
            assignin('base','cst',cst);
            this.handles=handles;
            changeWorkspace(this);
            
            % would need to remove
            %generateCstTable(this,cst);
            
        end
        function btObjRemove_Callback(this,hObject, ~)
            % hObject    handle to btnuiTableAdd (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles=this.handles;
            ix = get(hObject,'UserData');
            
            cst = evalin('base','cst');
            %Add Standard Objective
            
            cst{ix(1),6}(ix(2)) = [];
            
            assignin('base','cst',cst);
            this.handles=handles;
            changeWorkspace(this);
            
            %generateCstTable(this,cst);

        end
        
        function changeObjFunction_Callback(this,hObject, ~)
            % hObject    handle to btnuiTableAdd (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles=this.handles;
            data = get(hObject,'UserData');
            ix = data{1};
            classNames = data{2};
            classToCreate = classNames{get(hObject,'Value')};
            
            cst = evalin('base','cst');
            %Add Standard Objective
            
            %We just check if the user really wanted to change the objective to be
            %user-friendly
            currentObj = cst{ix(1),6}{ix(2)};
            currentClass = class(currentObj);
            if ~strcmp(currentClass,classToCreate)
                newObj = eval(classToCreate);
                
                % Only if we have a penalty value for optimization, apply the new one
                % Maybe this check should be more exact?
                if (isfield(currentObj,'penalty') || isprop(currentObj,'penalty')) && isprop(newObj,'penalty')
                    newObj.penalty = currentObj.penalty;
                end
                
                cst{ix(1),6}{ix(2)} = struct(newObj);
                
                assignin('base','cst',cst);
                this.handles=handles;
                changeWorkspace(this);
                
                %generateCstTable(this,cst);
            end
        end
        
        function editCstParams_Callback(this,hObject,~)
            handles=this.handles;
            data = hObject.UserData;
            ix = data(1);
            col = data(2);
            
            cst = evalin('base','cst');
            
            str = get(hObject,'String');
            val  = get(hObject,'Value');
            
            switch col
                case 2
                    cst{ix,col} = str;
                case 3
                    cst{ix,col} = str{val};
                case 5
                    cst{ix,col}.Priority = uint32(str2double(str));
                otherwise
                    warning('Wrong column assignment in GUI based cst setting');
            end
            
            assignin('base','cst',cst);
            this.handles=handles;
            changeWorkspace(this);
            
            %generateCstTable(this,cst);
        end
        
        function editObjParam_Callback(this,hObject, ~)
            % hObject    handle to btnuiTableAdd (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles=this.handles;
            ix = get(hObject,'UserData');
            
            cst = evalin('base','cst');
            
            if ix(3) == 0
                cst{ix(1),6}{ix(2)}.penalty = str2double(hObject.String);
            elseif isequal(hObject.Style,'popupmenu')
                cst{ix(1),6}{ix(2)}.parameters{ix(3)} = hObject.Value;
            else
                cst{ix(1),6}{ix(2)}.parameters{ix(3)} = str2double(hObject.String);
            end
            
            assignin('base','cst',cst);
            this.handles=handles;
            changeWorkspace(this);
            
            %generateCstTable(this,cst);
            
        end
    end
end
