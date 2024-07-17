classdef matRad_ParetoSliceWidget < matRad_Widget
    % matRad_ParetoSliceWidget class
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
        DosePlotAxes
    end        
    
    methods
        function this = matRad_ParetoSliceWidget(handleParent)
            this = this@matRad_Widget(handleParent);
            this.DosePlotAxes = axes(this.widgetHandle,'Position',[0 0 1 1]);
        end
        
        function this = initialize(this)
        end
                
        
    end
    
    methods (Access = protected)
        function this = createLayout(this)
                       
            parent = this.widgetHandle;
            
            matRad_cfg = MatRad_Config.instance();
                   

            this.createHandles();
        end
        
    end
    
end


