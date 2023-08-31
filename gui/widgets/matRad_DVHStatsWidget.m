classdef matRad_DVHStatsWidget < matRad_Widget

    % matRad_DVHStatsWidget class to generate GUI widget display DVH and
    % stats
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
    properties
        selectedDisplayOption;
        lockUpdate = false;
        dvhWidgetHandle = [];
        statWidgetHandle = [];
    end

    methods
        function this = matRad_DVHStatsWidget(selectedDisplayOption,handleParent)  

            
            matRad_cfg = MatRad_Config.instance();
            if nargin < 2 

                handleParent = figure(...
                    'Units','normalized',...
                    'OuterPosition',[0 0 0.5 1],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,...  'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','figure',...
                    'ToolBar','figure',...
                    'Name','MatRad Plan Analysis',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figDVHStat');
                
            end
            this = this@matRad_Widget(handleParent);
            this.selectedDisplayOption = selectedDisplayOption;
            this.dvhWidgetHandle.selectedCube = selectedDisplayOption;
            this.statWidgetHandle.selectedCube = selectedDisplayOption;
            this.lockUpdate = true;
            this.dvhWidgetHandle.lockUpdate = true;
            this.statWidgetHandle.lockUpdate = true;
            this.update();

        end        
        
        function this=update(this,evt)
            if this.lockUpdate
                if nargin == 2
                    doUpdate = this.checkUpdateNecessary({'resultGUI','cst','pln'},evt);

                    if doUpdate
                        this.dvhWidgetHandle.update(evt);
                        this.statWidgetHandle.update(evt);
                    end
                else
                    %Check for change in the selected cube 
                    if ~strcmp(this.dvhWidgetHandle.selectedCube, this.selectedDisplayOption)
                        this.dvhWidgetHandle.selectedCube = this.selectedDisplayOption;
                        this.statWidgetHandle.selectedCube = this.selectedDisplayOption;

                    end
                        this.dvhWidgetHandle = this.dvhWidgetHandle.update();
                        this.statWidgetHandle = this.statWidgetHandle.update();
                        %Clear previous DVH and stat
                        if numel(this.dvhWidgetHandle.widgetHandle.Children) > 2 
                            this.removeOverlap();
                        end
                 
                end
            end
        end
        function set.selectedDisplayOption(this,value)
            this.selectedDisplayOption = value;
            this.update();

        end
        function removeOverlap(this)
            % Clear previous plotted objects
            delete(this.dvhWidgetHandle.widgetHandle.Children(3)); 

        end
    end

    methods(Access = protected)
        function this = createLayout(this)
            h88 = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
            %DVH Panel
            p1 = uipanel(...
                'Parent',h88,...                      
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','panelDVH',...
                'Clipping','off',...
                'Position',[0.005 0.505 0.99 0.495],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Title','DVH');
            % Statistics panel
            p2 = uipanel(...
                'Parent',h88,...                
                'BackgroundColor',matRad_cfg.gui.backgroundColor,...
                'Tag','panelStats',...
                'Clipping','off',...
                'Position',[0.005 0.005 0.99 0.495],...
                'FontName',matRad_cfg.gui.fontName,...
                'FontSize',matRad_cfg.gui.fontSize,...
                'FontWeight',matRad_cfg.gui.fontWeight,...
                'Title','Statistics');
            
            %Initiate DVH and Stats Widgets
            this.dvhWidgetHandle = matRad_DVHWidget([],p1);
            this.statWidgetHandle = matRad_StatisticsWidget([],p2);
            this.createHandles();
            
        end
        
    end
end

