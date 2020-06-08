classdef matRad_StructureVisibilityWidget < matRad_Widget
    
    properties
        
    end
    
    methods
        function this = matRad_StructureVisibilityWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[250.4 45 30 15],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...  
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Structure Visibility',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
                
            end
            this = this@matRad_Widget(handleParent);
        end
        
        function this=initialize(this)
             this.update();
            
        end
        
        function this=update(this)
            if evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')
                updateStructureTable(this, evalin('base','cst'));
            else
                set(this.handles.legendTable,'String','');
            end
           
        end
        
        function changeWorkspace(obj)
           [env, ~] = matRad_getEnvironment();
            % handle environment
            switch env
                case 'MATLAB'
                    notify(obj, 'workspaceChanged');
                case 'OCTAVE'
            end
        end
    end
    
    methods (Access = protected)
        function this = createLayout(this)
            h86 = this.widgetHandle;
            
            h87 = uicontrol(...
                'Parent',h86,...
                'Units','normalized',...
                'Style','listbox',...
                'Value',1,...
                'Position',[0.02 0.01 0.97 0.98],...
                'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
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
            
            handles = this.handles;
            cst = evalin('base','cst');
            
            idx    = get(hObject,'Value');
            clr    = dec2hex(round(cst{idx,5}.visibleColor(:)*255),2)';
            clr    = ['#';clr(:)]';
            
            %Get the string entries
            tmpString = get(handles.legendTable,'String');
            
            if handles.VOIPlotFlag(idx)
                handles.VOIPlotFlag(idx) = false;
                cst{idx,5}.Visible = false;
                tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
            elseif ~handles.VOIPlotFlag(idx)
                handles.VOIPlotFlag(idx) = true;
                cst{idx,5}.Visible = true;
                tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"><center>&#10004;</center></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
            end
            set(handles.legendTable,'String',tmpString);
            
            % update cst in workspace accordingly
            assignin('base','cst',cst)
            this.handles = handles;
            changeWorkspace(this);
            %UpdatePlot(handles)
        end
        
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
                clr = dec2hex(round(cst{s,5}.visibleColor(:)*255),2)';
                clr = ['#';clr(:)]';
                if handles.VOIPlotFlag(s)
                    tmpString{s} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"><center>&#10004;</center></TD><TD>',cst{s,2},'</TD></TR> </table></html>'];
                else
                    tmpString{s} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"></TD><TD>',cst{s,2},'</TD></TR> </table></html>'];
                end
            end
            set(handles.legendTable,'String',tmpString);
            this.handles = handles;
        end
    end
end

