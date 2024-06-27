classdef matRad_StructureVisibilityWidget < matRad_Widget
    % matRad_StructureVisibilityWidget class to generate GUI widget to set
    % visibility of structure in viewing widget
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

    methods
        function this = matRad_StructureVisibilityWidget(handleParent)
            if nargin < 1
                matRad_cfg = MatRad_Config.instance();
                handleParent = figure(...
                      'Units','characters',...
                      'Position',[200 10 30 30],...
                      'Visible','on',...
                      'Color',matRad_cfg.gui.backgroundColor,...
                      'IntegerHandle','off',...
                      'MenuBar','none',...
                      'Name','MatRad Structure Visibility',...
                      'NumberTitle','off',...
                      'HandleVisibility','callback',...
                      'Tag','figure1');
            end
            this = this@matRad_Widget(handleParent);
        end
    end
   
    
    methods (Access = protected)
        function this = createLayout(this)
            h86 = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
            % List box of stuctures that can be selected for display

            if this.isInUifigure
                pos = getpixelposition(h86);
                h87 = uilistbox(h86);
                h87.Position = [0 0 pos(3) pos(4)];
                h87.Tooltip = 'Choose which structures should be displayed in the GUI';
                h87.BackgroundColor = matRad_cfg.gui.elementColor;
                h87.FontColor = matRad_cfg.gui.textColor;
                h87.FontSize = matRad_cfg.gui.fontSize;
                h87.FontWeight = matRad_cfg.gui.fontWeight;
                h87.FontName = matRad_cfg.gui.fontName;
                h87.ClickedFcn = @(hObject,eventdata) legendTable_Callback(this,hObject,eventdata);
                h87.Items = {'no data loaded'};
                h87.Tag = 'legendTable';
            else
                h87 = uicontrol(...
                    'Parent',h86,...
                    'Units','normalized',...
                    'TooltipString','Choose which structures should be displayed in the GUI',...
                    'HorizontalAlignment','left',...
                    'Max',1,...
                    'Style','listbox',...
                    'Value',1,...
                    'String','no data loaded',...
                    'Position',[0.02 0.01 0.97 0.98],...
                    'BackgroundColor',matRad_cfg.gui.elementColor,...
                    'ForegroundColor',matRad_cfg.gui.textColor,...
                    'FontSize',matRad_cfg.gui.fontSize,...
                    'FontWeight',matRad_cfg.gui.fontWeight,...
                    'FontName',matRad_cfg.gui.fontName,...
                    'Callback',@(hObject,eventdata) legendTable_Callback(this,hObject,eventdata),...
                    'Tag','legendTable');
            end
            this.createHandles();
         
        end

        function this=doUpdate(this,evt)
            if evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')
                updateStructureTable(this, evalin('base','cst'));
            else
                if this.isInUifigure()
                    this.handles.legendTable.Items = {'no data loaded'};
                else
                    set(this.handles.legendTable,'String','no data loaded');
                end
            end           
        end
    end

    
    
    methods (Access = protected)
        function legendTable_Callback(this, hObject, event)
            % hObject    handle to legendTable (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: contents = cellstr(get(hObject,'String')) returns legendTable contents as cell array
            %        contents{get(hObject,'Value')} returns selected item from legendTable
            
            if this.isInUifigure
                content = hObject.Items;
                idx = event.InteractionInformation.Item;
            else
                content = get(hObject,'String');
                idx = get(hObject,'Value');
            end
            if numel(content) == 1 && strcmp(content,'no data loaded')
                return;
            end
            
            handles = this.handles;
            cst = evalin('base','cst');
            
            if handles.VOIPlotFlag(idx)
                handles.VOIPlotFlag(idx) = false;
                cst{idx,5}.Visible = false;
            else
                handles.VOIPlotFlag(idx) = true;
                cst{idx,5}.Visible = true;
            end
              
            this.updateStructureTable(cst);
            
            % update cst in workspace accordingly
            assignin('base','cst',cst)
            this.handles = handles;
            this.changedWorkspace('cst');
            %UpdatePlot(handles)
        end
        
        %Update cst with Visibility settings
        function cst = updateStructureTable(this,cst)
            handles=this.handles;
            colorAssigned = true;

            
            % check whether all structures have an assigned color
            for i = 1:size(cst,1)
                if ~isfield(cst{i,5},'visibleColor')
                    colorAssigned = false;
                    break;
                elseif isempty(cst{i,5}.visibleColor)
                    colorAssigned = false;
                    break;
                end
            end
            
            % assign color if color assignment is not already present or inconsistent
            if colorAssigned == false
                m         = 64;
                colorStep = ceil(m/size(cst,1));
                colors    = colorcube(colorStep*size(cst,1));
                % spread individual VOI colors in the colorcube color palette
                colors    = colors(1:colorStep:end,:);
                
                for i = 1:size(cst,1)
                    cst{i,5}.visibleColor = colors(i,:);
                end
            end
            
            for s = 1:size(cst,1)
                handles.VOIPlotFlag(s) = cst{s,5}.Visible;
                
                [tmpString{s},tmpStyles{s}] = this.getListEntry(cst(s,:));
            end
            if this.isInUifigure
                handles.legendTable.Items = tmpString;
                for s = 1:numel(tmpStyles)
                    addStyle(handles.legendTable,tmpStyles{s},'Item',s);
                end
            else
                set(handles.legendTable,'String',tmpString);
            end
            this.handles = handles;
        end

        function [item,style] = getListEntry(this,cstElement)
            clr = cstElement{1,5}.visibleColor;

            if cstElement{1,5}.Visible
                checkbox = '☑ ';
            else
                checkbox = '☐ ';
            end

            % html is not supported in octave
            if this.isInUifigure()
                % calculate text color
                intensity = clr * [0.299 0.587 0.114]';
                if intensity > 150/255 %186/255 %https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color
                    txtClr = [0 0 0];
                else
                    txtClr = [1 1 1];
                end

                item = [checkbox cstElement{1,2}];
                style = uistyle('BackgroundColor',clr,'FontColor',txtClr);
            else
                matRad_cfg = MatRad_Config.instance();
                switch matRad_cfg.env
                    case 'OCTAVE'
                        item = [checkbox cstElement{1,2}];
                    otherwise
                        hexClr = dec2hex(round(cstElement{1,5}.visibleColor(:)*255),2)';
                        hexClr = ['#';hexClr(:)]';
                        if cstElement{1,5}.Visible
                            item = ['<html><table border=0 ><TR><TD bgcolor=',hexClr,' width="18"><center>&#10004;</center></TD><TD>',cstElement{1,2},'</TD></TR> </table></html>'];
                        else
                            item = ['<html><table border=0 ><TR><TD bgcolor=',hexClr,' width="18"></TD><TD>',cstElement{1,2},'</TD></TR> </table></html>'];
                        end
                end
                style = [];
            end
        end
    end
end

