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
  
        function this=update(this,evt)
            if evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')
                updateStructureTable(this, evalin('base','cst'));
            else
                set(this.handles.legendTable,'String','no data loaded');
            end
           
        end
        
    end
   
    
    methods (Access = protected)
        function this = createLayout(this)
            h86 = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
            % List box of stuctures that can be selected for display
            h87 = uicontrol(...
                'Parent',h86,...
                'Units','normalized',...
                'Tooltip','Choose which structures should be displayed in the GUI',...
                'HorizontalAlignment','left',...
                'Max',1,...
                'Style','listbox',...
                'Value',1,...
                'Position',[0.02 0.01 0.97 0.98],...
                'BackgroundColor',matRad_cfg.gui.elementColor,...
                'ForegroundColor',matRad_cfg.gui.textColor,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'FontName',matRad_cfg.gui.fontName,...
                'Callback',@(hObject,eventdata) legendTable_Callback(this,hObject,eventdata),...
                'Tag','legendTable');
            
            this.createHandles();
         
        end
    end
    
    methods (Access = protected)
        function legendTable_Callback(this, hObject, event)
            % hObject    handle to legendTable (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            
            % Hints: contents = cellstr(get(hObject,'String')) returns legendTable contents as cell array
            %        contents{get(hObject,'Value')} returns selected item from legendTable
            
            if strcmp(get(hObject,'String'),'no data loaded')
                return;
            end
            
            handles = this.handles;
            cst = evalin('base','cst');
            
            idx    = get(hObject,'Value');
            clr    = dec2hex(round(cst{idx,5}.visibleColor(:)*255),2)';
            clr    = ['#';clr(:)]';
            
            %Get the string entries
            tmpString = get(handles.legendTable,'String');
            
            matRad_cfg = MatRad_Config.instance();
            
            % html not supported in octave       
            if handles.VOIPlotFlag(idx)
                handles.VOIPlotFlag(idx) = false;
                cst{idx,5}.Visible = false;
                switch matRad_cfg.env
                    case 'OCTAVE'
                        tmpString{idx} = ['☐ ' cst{idx,2}];
                    otherwise
                        tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
                end
            elseif ~handles.VOIPlotFlag(idx)
                handles.VOIPlotFlag(idx) = true;
                cst{idx,5}.Visible = true;
                switch matRad_cfg.env
                    case 'OCTAVE'
                        tmpString{idx} = ['☑ ' cst{idx,2}];
                    otherwise
                        tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"><center>&#10004;</center></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
                end
            end
           
            set(handles.legendTable,'String',tmpString);
            
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
            
            matRad_cfg = MatRad_Config.instance();
            
            for s = 1:size(cst,1)
                handles.VOIPlotFlag(s) = cst{s,5}.Visible;
                clr = dec2hex(round(cst{s,5}.visibleColor(:)*255),2)';
                clr = ['#';clr(:)]';
                % html is not supported in octave 
                
                switch matRad_cfg.env
                    case 'OCTAVE'
                        if handles.VOIPlotFlag(s)
                            tmpString{s} = ["o " cst{s,2}];
                        else
                            tmpString{s} = ["x " cst{s,2}];
                        end
                    otherwise
                        if handles.VOIPlotFlag(s)
                            tmpString{s} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"><center>&#10004;</center></TD><TD>',cst{s,2},'</TD></TR> </table></html>'];
                        else
                            tmpString{s} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"></TD><TD>',cst{s,2},'</TD></TR> </table></html>'];
                        end
                end
                
                
            end
            set(handles.legendTable,'String',tmpString);
            this.handles = handles;
        end
    end
end

