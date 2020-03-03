classdef matRad_OptimizationWidget < matRad_Widget
    
    properties
        
    end
    
    methods
        function this = matRad_OptimizationWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[170.4 45 150.4 35.5384615384615],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Optimization',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
                
            end
            this = this@matRad_Widget(handleParent);
        end
        
    end
    
    methods (Access = protected)
        function this = createLayout(this)
            h1 = this.widgetHandle;
            
            
        end
    end
    
    methods
        function cst = generateCstTable(handles,cst)
            handles = this.handles;
            cst = updateStructureTable(handles,cst);
            cstPanel = handles.uipanel3;
            
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
            
            %columnname = {'VOI name','VOI type','priority','obj. / const.'};%,'penalty','dose', 'EUD','volume','robustness'};
            
            %Get all Classes & classNames
            classNames = matRad_getObjectivesAndConstraints();
            % Collect Class-File & Display Names
            %classNames = {classList.Name; p.DefaultValue};
            
            %columnformat = {cst(:,2)',{'OAR','TARGET'},'numeric',...
            %       AllObjectiveFunction,...
            %       'numeric','numeric','numeric','numeric',{'none','WC','prob'}};
            
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
                        
                        %VOI
                        %data{Counter,1}  = cst{i,2};
                        %ypos = cstPanelPos(4) - (yTopSep + cnt*lineHeight);
                        xPos = 0.01;%5;
                        
                        %h = uicontrol(cstPanel,'Style','popupmenu','String',cst(:,2)','Position',[xPos ypos 100 objHeight]);
                        %h.Value = i;
                        h = uicontrol(cstPanel,'Style','pushbutton','String','-','Units','normalized','Position',[xPos ypos(cnt) buttonW objHeight],'TooltipString','Remove Objective/Constraint','Callback',{@btObjRemove_Callback,handles},...
                            'UserData',[i,j]);
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        h = uicontrol(cstPanel','Style','edit','String',cst{i,2},'Units','normalized','Position',[xPos ypos(cnt) nameW objHeight],'TooltipString','Name',...
                            'Enable','inactive',... %Disable editing of name atm
                            'UserData',[i,2],'Callback',{@editCstParams_Callback,handles}); %Callback added, however, editing is disabled atm
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        h = uicontrol(cstPanel,'Style','popupmenu','String',organTypes','Value',find(strcmp(cst{i,3},organTypes)),'Units','normalized','Position',[xPos ypos(cnt) typeW objHeight],'TooltipString','Segmentation Classification',...
                            'UserData',[i,3],'Callback',{@editCstParams_Callback,handles});
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        h = uicontrol(cstPanel,'Style','edit','String',num2str(cst{i,5}.Priority),'Units','normalized','Position',[xPos ypos(cnt) opW objHeight],'TooltipString',['Overlap Priority' newline '(Smaller number overlaps higher number)'],...
                            'UserData',[i,5],'Callback',{@editCstParams_Callback,handles});
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        
                        h = uicontrol(cstPanel,'Style','popupmenu','String',classNames(2,:)','Value',find(strcmp(obj.name,classNames(2,:))),'Units','normalized','Position',[xPos ypos(cnt) functionW objHeight],'TooltipString','Select Objective/Constraint',...
                            'UserData',{[i,j],classNames(1,:)},'Callback',{@changeObjFunction_Callback,handles});
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        
                        %Check if we have an objective to display penalty
                        if isa(obj,'DoseObjectives.matRad_DoseObjective')
                            h = uicontrol(cstPanel,'Style','edit','String',num2str(obj.penalty),'Units','normalized','Position',[xPos ypos(cnt) penaltyW objHeight],'TooltipString','Objective Penalty','UserData',[i,j,0],'Callback',{@editObjParam_Callback,handles});
                        else
                            h = uicontrol(cstPanel,'Style','edit','String','----','Units','normalized','Position',[xPos ypos(cnt) penaltyW objHeight],'Enable','off');
                        end
                        tmp_pos = get(h,'Position');
                        xPos = xPos + tmp_pos(3) + fieldSep;
                        
                        for p = 1:numel(obj.parameterNames)
                            %h = uicontrol(cstPanel,'Style','edit','String',obj.parameters{1,p},'Position',[xPos ypos(cnt) 100 objHeight],'Enable','inactive');
                            %xPos = xPos + h.Position(3) + fieldSep;
                            h = text('Parent',tmpAxes,'String',['| ' obj.parameterNames{p} ':'],'VerticalAlignment','middle','Units','normalized','Position',[xPos ypos(cnt)+lineHeight/2],'Interpreter','tex','FontWeight','normal',...
                                'FontSize',get(cstPanel,'FontSize'),'FontName',get(cstPanel,'FontName'),'FontUnits',get(cstPanel,'FontUnits'),'FontWeight','normal');%[xPos ypos(cnt) 100 objHeight]);
                            tmp_pos = get(h,'Extent');
                            xPos = xPos + tmp_pos(3) + fieldSep;
                            %h = annotation(cstPanel,'textbox','String',obj.parameters{1,p},'Units','pix','Position', [xPos ypos(cnt) 100 objHeight],'Interpreter','Tex');
                            
                            %Check if we have a cell and therefore a parameter list
                            if iscell(obj.parameterTypes{p})
                                h = uicontrol(cstPanel,'Style','popupmenu','String',obj.parameterTypes{p}','Value',obj.parameters{p},'TooltipString',obj.parameterNames{p},'Units','normalized','Position',[xPos ypos(cnt) paramW*2 objHeight],'UserData',[i,j,p],'Callback',{@editObjParam_Callback,handles});
                            else
                                h = uicontrol(cstPanel,'Style','edit','String',num2str(obj.parameters{p}),'TooltipString',obj.parameterNames{p},'Units','normalized','Position',[xPos ypos(cnt) paramW objHeight],'UserData',[i,j,p],'Callback',{@editObjParam_Callback,handles});
                            end
                            
                            tmp_pos = get(h,'Extent');
                            xPos = xPos + tmp_pos(3) + fieldSep;
                        end
                        
                        cnt = cnt +1;
                    end
                end
            end
            xPos = 0.01; %5
            hAdd = uicontrol(cstPanel,'Style','pushbutton','String','+','Units','normalized','Position',[xPos ypos(cnt) buttonW objHeight],'TooltipString','Add Objective/Constraint','Callback',{@btObjAdd_Callback,handles});
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
        
    end
end
