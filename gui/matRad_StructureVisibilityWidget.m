classdef matRad_StructureVisibilityWidget < matRad_Widget
    
    properties
        
    end
    
    methods
        function this = matRad_StructureVisibilityWidget(handleParent)
            if nargin < 1
                handleParent = figure(...
                    'Units','characters',...
                    'Position',[138.4 -7.38461538461539 273.4 59.5384615384615],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...   'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','matRadGUI',...
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
            
            % guidata(hObject, handles);
            this.handles = handles;
            UpdatePlot(handles)
        end
    end
end

