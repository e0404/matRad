classdef matRad_DVHStatsWidget < matRad_Widget
    
    properties
        SelectedCube;
        
        dvhWidgetHandle = [];
        statWidgetHandle = [];
    end
    
    events
        
    end
    
    methods
        function this = matRad_DVHStatsWidget(varargin)    
            
            p = inputParser;
            p.addOptional('handleParent',[], @ishandle);
            p.addParameter('SelectedCube','physicalDose', @ischar);
            p.parse(varargin{:});
            
            matRad_cfg = MatRad_Config.instance();
            if isempty(p.Results.handleParent)
                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.005 0.05 0.495 0.9],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,...  'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','figure',...
                    'ToolBar','figure',...
                    'Name','MatRad Plan Analysis',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1');
                
            end
            
            this = this@matRad_Widget(handleParent);   
            this.SelectedCube = p.Results.SelectedCube;    
            this.dvhWidgetHandle.SelectedCube = this.SelectedCube;
            this.statWidgetHandle.SelectedCube = this.SelectedCube;
            this.showDVHstat();
        end        
        
        function this=update(this,evt)
            if nargin == 2
                doUpdate = this.checkUpdateNecessary({'resultGUI','cst','pln'},evt);
                
                if doUpdate
                    this.dvhWidgetHandle.update(evt);
                    this.statWidgetHandle.update(evt);
                end
            else                  
                this.dvhWidgetHandle.update();
                this.statWidgetHandle.update();
            end
        end
        
    end
    
    methods(Access = protected)
        function this = createLayout(this)
            h88 = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
            
            p1 = uipanel(...
                'Parent',h88,...                      
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel1',...
                'Clipping','off',...
                'Position',[0.005 0.505 0.99 0.495],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Title','DVH');
            
            p2 = uipanel(...
                'Parent',h88,...                
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','uipanel2',...
                'Clipping','off',...
                'Position',[0.005 0.005 0.99 0.495],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Title','Statistics');
            
            this.dvhWidgetHandle = matRad_DVHWidget(p1);
            this.statWidgetHandle = matRad_StatisticsWidget(p2);
            
            this.createHandles();
            
        end
        
        function this = showDVHstat(this)
            this.dvhWidgetHandle.showDVH();
            this.statWidgetHandle.showStatistics();
        end 
    end
end