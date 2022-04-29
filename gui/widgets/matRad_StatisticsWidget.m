classdef matRad_StatisticsWidget < matRad_Widget
    
    properties
        SelectedCube='physicalDose';
    end
    
    events
        
    end
    
    methods
        function this = matRad_StatisticsWidget(varargin)
            
            p = inputParser;
            p.addOptional('handleParent',[], @ishandle);
            p.addParameter('SelectedCube','physicalDose', @ischar);
            p.parse(varargin{:});
            
            if isempty(p.Results.handleParent)
                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.005 0.05 0.495 0.45],...
                    'Visible','on',...
                    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...  'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Statistics',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1',...
                    'PaperSize',[20.99999864 29.69999902]);
            else 
                handleParent  = p.Results.handleParent;               
                
            end
            
            this = this@matRad_Widget(handleParent);   
            this.SelectedCube = p.Results.SelectedCube;

        end

        
        function this=update(this)
            doUpdate = true;
            if nargin == 2
                doUpdate = this.checkUpdateNecessary({'resultGUI','cst','pln'},evt);
            end
            
            if ~doUpdate
                return;
            end
            
            if evalin('base','exist(''resultGUI'')')
                this.showStatistics();
            end
        end
        
    end
    
    methods(Access = protected)
        function this = createLayout(this)
            h88 = this.widgetHandle;           
            this.createHandles();
            
        end
    end
    
    methods
         function showStatistics(this)
            
            resultGUI = evalin('base','resultGUI');
            %Content = get(handles.popupDisplayOption,'String');
            
            %SelectedCube = Content{get(handles.popupDisplayOption,'Value')};
            
            pln = evalin('base','pln');
            resultGUI_SelectedCube.physicalDose = resultGUI.(this.SelectedCube);
            
            if ~strcmp(pln.propOpt.bioOptimization,'none')
                
                %check if one of the default fields is selected
                if sum(strcmp(this.SelectedCube,{'physicalDose','effect','RBE,','RBExDose','alpha','beta'})) > 0
                    resultGUI_SelectedCube.physicalDose = resultGUI.physicalDose;
                    resultGUI_SelectedCube.RBExDose     = resultGUI.RBExDose;
                else
                    Idx    = find(this.SelectedCube == '_');
                    SelectedSuffix = this.SelectedCube(Idx(1):end);
                    resultGUI_SelectedCube.physicalDose = resultGUI.(['physicalDose' SelectedSuffix]);
                    resultGUI_SelectedCube.RBExDose     = resultGUI.(['RBExDose' SelectedSuffix]);
                end
            end
            
            cst = evalin('base','cst');
            %matRad_indicatorWrapper(this.widgetHandle,cst,pln,resultGUI_SelectedCube);
            
            if isfield(resultGUI_SelectedCube,'RBExDose')
                doseCube = resultGUI_SelectedCube.RBExDose;
            else
                doseCube = resultGUI_SelectedCube.physicalDose;
            end
            
            if ~exist('refVol', 'var')
                refVol = [];
            end
            
            if ~exist('refGy', 'var')
                refGy = [];
            end
            
            qi  = matRad_calcQualityIndicators(cst,pln,doseCube,refGy,refVol);

            ixVoi = cellfun(@(c) c.Visible == 1,cst(:,5));
            qi = qi(ixVoi);
            matRad_showQualityIndicators(this.widgetHandle,qi);

         end
        
    end
end