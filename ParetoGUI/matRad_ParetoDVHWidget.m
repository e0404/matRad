classdef matRad_ParetoDVHWidget < matRad_Widget
    % matRad_ParetoDVHWidget class creates a widget for the Pareto GUI to
    % display plan DVH's
    % 
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024 the matRad development team. 
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
        DVHPlotAxes
    end        
    
    methods
        function this = matRad_ParetoDVHWidget(handleParent)
            this = this@matRad_Widget(handleParent);
            this.DVHPlotAxes = axes(this.widgetHandle,'Position',[0 0 1 1]);
            this.plotInitialDVH();
        end
        
        function this = initialize(this)
            this.update();
        end
    
    end
    
    methods (Access = protected)
        function this = createLayout(this)
                       
            parent = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
                   

            this.createHandles();
        end
        
        function this = doUpdate(this,evt)
            %getFromWorkspace(this);
            %updateInWorkspace(this);
        end
    end

    methods (Access = private)
        function plotInitialDVH(this)
            ParetoHelperObject = evalin('base','ParetoHelperObject');
            %Shows the DVH for the current plan
            resultGUI = matRad_calcCubes(ParetoHelperObject.currentWeights,evalin('base','dij'));
            dvh = matRad_calcDVH(evalin('base','cst'),resultGUI.physicalDose,'cum');
    
            matRad_showDVH(this.DVHPlotAxes,...
                dvh,evalin('base','cst'),evalin('base','pln'),ParetoHelperObject.linestyle);
            if ParetoHelperObject.linestyle < 4
                ParetoHelperObject.linestyle = ParetoHelperObject.linestyle + 1;
            else
                ParetoHelperObject.linestyle = 1;
            end
        end
    end
    
end


