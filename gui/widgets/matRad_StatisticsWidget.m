classdef matRad_StatisticsWidget < matRad_Widget
    % matRad_StatisticsWidget class to generate GUI widget to display plan
    % statistics.
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
        selectedCube;
        lockUpdate = false;

    end
    
    events
        
    end
    
    methods

        function this = matRad_StatisticsWidget(SelectedCube,handleParent)    % use (varargin) ?
             
             matRad_cfg = MatRad_Config.instance();
             if nargin < 2 

                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.005 0.05 0.495 0.45],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,...  'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad Statistics',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figStat',...
                    'PaperSize',[20.99999864 29.69999902]);
        
                
              end
              this = this@matRad_Widget(handleParent);
              this.selectedCube = SelectedCube;
        end

        
        function this=update(this)

            if this.lockUpdate
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
    end
    
    methods(Access = protected)
        function this = createLayout(this)
            h88 = this.widgetHandle;           
            this.createHandles();
            
        end
    end
    
    methods
         function set.selectedCube(this,value)
            this.selectedCube=value;
        end
         function showStatistics(this)
            
            resultGUI = evalin('base','resultGUI');
            pln = evalin('base','pln');
            cst = evalin('base','cst');
            doseCube = resultGUI.(this.selectedCube);
            
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