classdef matRad_DVHWidget < matRad_Widget
    % matRad_DVHWidget class to generate GUI widget to display DVH 
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
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        selectedCube = [];
    end

    properties (SetAccess = private)
        dvhAx;
    end
      
    methods
        function this = matRad_DVHWidget(handleParent) 
           
            matRad_cfg = MatRad_Config.instance();
            if nargin < 1 
                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.005 0.5 0.495 0.45],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,...  'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad DVH',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figDVH',...
                    'PaperSize',[20.99999864 29.69999902]);
                           
            end
            this = this@matRad_Widget(handleParent);             
        end
                       
        function removeOverlap(this)
            %Clear previous plotted objects from the figure
            delete(this.widgetHandle.Children(3));
        end

        function initialize(this)
            initialize@matRad_Widget(this);
        end
    end
        
    methods(Access = protected)
        function this = createLayout(this,handleParent)
            this.createHandles();            
        end

        function this = doUpdate(this,evt)
            if ~this.updateLock && ~isempty(this.selectedCube)
                doUpdate = true;
                if nargin == 2
                    doUpdate = this.checkUpdateNecessary({'resultGUI','cst','pln'},evt);
                end

                if doUpdate && evalin('base','exist(''resultGUI'')') && evalin('base','exist(''cst'')')
                    this.showDVH();
                    ch = get(this.widgetHandle,'Children'); %Octave compatibility, no dot indexing
                    if numel(ch) > 2
                        this.removeOverlap();
                    end
                end
            end
        end
    end
    
    methods

        function set.selectedCube(this,value)
            this.selectedCube=value;
            this.update();
        end

        function showDVH(this)
            if isempty(this.selectedCube)
                return;
            end

            matRad_cfg = MatRad_Config.instance();
            
            %Create axis if not existing, otherwise clear
            if isgraphics(this.dvhAx)
                cla(this.dvhAx);
            else
                this.dvhAx = axes(this.widgetHandle,...
                    'Color',matRad_cfg.gui.elementColor,...
                    'XColor',matRad_cfg.gui.textColor,...
                    'YColor',matRad_cfg.gui.textColor,...
                    'GridColor',matRad_cfg.gui.textColor,...
                    'MinorGridColor',matRad_cfg.gui.backgroundColor);
            end
            resultGUI = evalin('base','resultGUI');
            pln = evalin('base','pln');
            cst = evalin('base','cst');
            ct  = evalin('base','ct');

            % Calculate and show DVH
            doseCube = resultGUI.(this.selectedCube);
            dvh = matRad_calcDVH(cst,doseCube,'cum');
            
            matRad_showDVH(dvh,cst,pln,'axesHandle',this.dvhAx);

            if ~isfield(pln,'multScen')
                multScen = matRad_NominalScenario(ct);
            else
                multScen = pln.multScen;
            end
            multScen = matRad_ScenarioModel.create(multScen,ct);
            
            %check scenarios
            if multScen.totNumScen > 1
                for i = 1:multScen.totNumScen
                    scenFieldName = sprintf('%s_scen%d',this.selectedCube,i);
                    if isfield(resultGUI,scenFieldName)
                        tmpDvh = matRad_calcDVH(cst,resultGUI.(scenFieldName),'cum'); % Calculate cumulative scenario DVH
                        matRad_showDVH(tmpDvh,cst,pln,'axesHandle',this.dvhAx,'plotLegend',false,'LineWidth',0.5,'LineStyle','--'); % Show DVH plot
                    end
                end
            end
            
            %No dot indexing for Octave compatibility
            hTitle = get(this.dvhAx,'Title');
            set(hTitle,'String', strrep(this.selectedCube, '_',' '), 'Color',matRad_cfg.gui.highlightColor);
        end
        
    end
end