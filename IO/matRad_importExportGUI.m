classdef matRad_importExportGUI
   
    properties
        guiHandle
    end
    
    properties (Access = private)
        importWidget
        exportWidget
    end
        
    
    methods
        function obj = matRad_importExportGUI()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.guiHandle = figure;
            p1 = uipanel(obj.guiHandle,'Position',[0 0 0.5 1]);
            p2 = uipanel(obj.guiHandle,'Position',[0.5 0 0.5 1]);
            
            obj.importWidget = matRad_importWidget(p1);
            obj.exportWidget = matRad_exportWidget(p2);
        end
        
    end
end

